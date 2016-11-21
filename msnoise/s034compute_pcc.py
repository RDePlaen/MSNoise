""" This code is responsible for the computation of the cross-correlation
functions.

This script will group *jobs* marked "T"odo in the database by day and process
them using the following scheme. As soon as one day is selected, the
corresponding jobs are marked "I"n Progress in the database. This allows
running several instances of this script in parallel.

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* |cc_sampling_rate|
* |analysis_duration|
* |overlap|
* |maxlag|
* |corr_duration|
* |windsorizing|
* |resampling_method|
* |decimation_factor|
* |remove_response|
* |response_format|
* |response_path|
* |response_prefilt|
* |preprocess_lowpass|
* |preprocess_highpass|
* |keep_all|
* |keep_days|
* |stack_method| | *new in 1.4*
* |pws_timegate| | *new in 1.4*
* |pws_power| | *new in 1.4*


Waveform Pre-processing
~~~~~~~~~~~~~~~~~~~~~~~
Pairs are first split and a station list is created. The database is then
queried to get file paths. For each station, all files potentially containing
data for the day are opened. The traces are then merged and splitted, to obtain
the most continuous chunks possible. The different chunks are then demeaned,
tapered and merged again to a 1-day long trace. If a chunk is not aligned
on the sampling grid (that is, start at a integer times the sample spacing in s)
, the chunk is phase-shifted in the frequency domain. This requires tapering and
fft/ifft. If the gap between two chunks is small, compared to a currently
hard-coded value (10 samples), the gap is filled with interpolated values.
Larger gaps will not be filled with interpolated values and remaining chunks
will be tapered and then merged with 0 values in the gaps.

If shorter than 1-day, the trace final is padded with zeros. If longer, it is
cut to match the start/end of the day.

If configured, each 1-day long trace is corrected for its instrument response.
Currently, only dataless seed and inventory XML are supported.

.. note:: Removing the instrument response is a computationally very expensive
   task and *not* useful for dv/v iff your instruments didn't change during the
   analysed period. It is also not needed for tomography iff all instruments are
   the same, or at least have an identical phase response in the frequency band
   of interest.

Each 1-day long trace is then low-passed (at ``preprocess_lowpass`` Hz),
high-passed (at ``preprocess_highpass`` Hz), then if needed,
decimated/downsampled. Decimation/Downsampling are configurable
(``resampling_method``) and users are advised testing both. One advantage of
Downsampling over Decimation is that it is able to downsample the data by any
factor, not only integer factors.

.. note:: Python 3 users will most probably struggle installing
    scikits.samplerate, and therefore will have to use either Decimate or
    Lanczos instead of Resample. This is not a problem because the Lanczos
    resampling give very similar results as the scikits.samplerate package.

Processing
~~~~~~~~~~

Once all traces are preprocessed, station pairs are processed sequentially.
If a component different from *ZZ* is to be computed, the traces are first
rotated. This supposes the user has provided the station coordinates in the
*station* table. The rotation is computed for Radial and Transverse components:

.. code-block:: python

    R = N * np.cos(Az * np.pi / 180.) + E * np.sin(Az * np.pi / 180.)
    T = N * np.sin(Az * np.pi / 180.) - E * np.cos(Az * np.pi / 180.)

Then, for each ``corr_duration`` window in the signal, and for each filter
configured in the database, the traces are clipped to ``windsorizing`` times
the RMS (or 1-bit converted) and then whitened (see :ref:`whiten`) between the
frequency bounds.
When both traces are ready, the cross-correlation function is computed
(see :ref:`mycorr`). The function returned contains data for time lags
corresponding to ``maxlag`` in the acausal (negative lags) and causal
(positive lags) parts.

Stacking and Saving Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If configured (setting ``keep_all`` to 'Y'), each ``corr_duration`` CCF is
saved to the hard disk. By default, the ``keep_days`` setting is set to True
and so "N = 1 day / corr_duration" CCF are stacked and saved to the hard disk
in the STACKS/001_DAYS folder.

.. note:: Currently, the keep-all data (every CCF) are not used by next steps.

If ``stack_method`` is 'linear', then a simple mean CFF of all windows is saved
as the daily CCF. On the other hand, if ``stack_method`` is 'pws', then
all the Phase Weighted Stack (PWS) is computed and saved as the daily CCF. The
PWS is done in two steps: first the mean coherence between the instataneous
phases of all windows is calculated, and eventually serves a weighting factor
on the mean. The smoothness of this weighting array is defined using the
``pws_timegate`` parameter in the configuration. The weighting array is the
power of the mean coherence array. If ``pws_power`` is equal to 0, a linear
stack is done (then it's faster to do set ``stack_method`` = 'linear'). Usual
value is 2.

.. warning:: PWS is largely untested, not cross-validated. It looks good, but
    that doesn't mean a lot, does it? Use with Caution! And if you
    cross-validate it, please let us know!!

.. seealso:: Schimmel, M. and Paulssen H., "Noise reduction and detection
    of weak, coherent signals through phase-weighted stacks". Geophysical Journal
    International 130, 2 (1997): 497-505.

Once done, each job is marked "D"one in the database.

To run this script:

.. code-block:: sh

    $ msnoise compute_ac


This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 compute_ac

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.


.. versionadded:: 1.4
    The Instrument Response removal & The Phase Weighted Stack &
    Parallel Processing

.. versionadded:: 1.5
    The Obspy Lanczos resampling method, gives similar results as the
    scikits.samplerate package, thus removing the requirement for it.

"""
import sys
import time
from scipy.fftpack.helper import next_fast_len

try:
    from scikits.samplerate import resample
except:
    pass

from .api import *
#from .move2obspy import myCorr
#from .move2obspy import whiten
from .Pnowhiten import nowhiten
from .myPCorr import myPCorr
from .preprocessing import preprocess


class Params():
    pass


def main():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('*** Starting: Compute AC ***')

    # Connection to the DB
    db = connect()

    if len(get_filters(db, all=False)) == 0:
        logging.info("NO FILTERS DEFINED, exiting")
        sys.exit()

    # Get Configuration
    params = Params()
    params.goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    params.goal_duration = float(get_config(db, "analysis_duration"))
    params.overlap = float(get_config(db, "overlap"))
    params.maxlag = float(get_config(db, "maxlag"))
    params.min30 = float(get_config(db, "corr_duration")) * params.goal_sampling_rate
    params.windsorizing = float(get_config(db, "windsorizing"))
    params.resampling_method = get_config(db, "resampling_method")
    params.decimation_factor = int(get_config(db, "decimation_factor"))
    params.preprocess_lowpass = float(get_config(db, "preprocess_lowpass"))
    params.preprocess_highpass = float(get_config(db, "preprocess_highpass"))
    params.keep_all = get_config(db, 'keep_all', isbool=True)
    params.keep_days = get_config(db, 'keep_days', isbool=True)
#Raph    params.components_to_compute = get_components_to_compute(db)

    params.stack_method = get_config(db, 'stack_method')
    params.pws_timegate = float(get_config(db, 'pws_timegate'))
    params.pws_power = float(get_config(db, 'pws_power'))

    params.components_to_compute = ['Z', 'E', 'N']
    ##################################

    logging.info("Will compute %s" % " ".join(params.components_to_compute))
    ##################################
    stations_to_analyse = ["%s.%s" % (sta.net, sta.sta) for sta in get_stations(db, all=True)]#extract all stations

    pairs = []#pair of components
    ##################################
    #modified part to make pair list with comps
    ##################################
    i = 0
    for comp in params.components_to_compute:
        for newcomp in params.components_to_compute:
            if comp == newcomp:
                if i == 0:
                    pairs = np.array(':'.join([comp, newcomp]))
                    i+=1
                else:
                    pairs = np.vstack((pairs,':'.join([comp, newcomp])))
    pairs = np.hstack(pairs)
    print('Pairs = ', len(pairs))
    while is_next_job(db, jobtype='AC'):
        jobs = get_next_job(db, jobtype='AC')
        stations = []
        #Raph pairs = []
        refs = []

#Raph changed to fit AC
#        for job in jobs:
#            refs.append(job.ref)
#            pairs.append(job.pair)
#            netsta1, netsta2 = job.pair.split(':')
#            stations.append(netsta1)
#            stations.append(netsta2)
#            goal_day = job.day

        #go through job to make job array with stations
        for job in jobs:
            refs.append(job.ref)
            #pairs.append(job.pair)#find a way to pair comps   /!\
            netsta = job.pair  #just 1 station in the cell?
            stations.append(netsta)
            #stations.append(netstacomp2)
            goal_day = job.day



        stations = np.unique(stations)

        logging.info("New AC Job: %s (%i stations with %i pairs each)" %
                     (goal_day , len(stations), len(pairs)))
        jt = time.time()

        xlen = int(params.goal_duration * params.goal_sampling_rate)

#Raph changed here the trame selection

        comps = ['Z', 'E', 'N']
        tramef_Z = np.zeros((len(stations), xlen))
        tramef_E = np.zeros((len(stations), xlen))
        tramef_N = np.zeros((len(stations), xlen))
        basetime, tramef_Z, tramef_E, tramef_N = preprocess(db, stations, comps, goal_day, params, tramef_Z, tramef_E, tramef_N)# preprocessing


        # print '##### STREAMS ARE ALL PREPARED AT goal Hz #####'
        dt = 1. / params.goal_sampling_rate

        begins = []
        ends = []
        i = 0
        while i <=  (params.goal_duration - params.min30/params.goal_sampling_rate):
            begins.append(int(i * params.goal_sampling_rate))
            ends.append(int(i * params.goal_sampling_rate + params.min30))
            i += int(params.min30/params.goal_sampling_rate * (1.0-params.overlap))
        # ##########################################################################################################

        for station in stations:
            orig_pair = station
            for pair in pairs:
        # ITERATING OVER PAIRS #####
        #for pair in pairs:
                #orig_pair = pair
    
                #logging.info('Processing pair: %s' % pair.replace(':', ' vs '))
                logging.info("Processing pair %s for station %s" %(pair, station))
                tt = time.time()
                comp1,comp2=pair.split(':')
                components=comp1+comp2
                ### load trames
                #assign trames according to pair

                station_to_analyse=np.where(stations==station)

                if pair.split(':')[0]=='Z':
                    tr1=tramef_Z[station_to_analyse]
                elif pair.split(':')[0]=='E':
                    tr1=tramef_E[station_to_analyse]
                elif pair.split(':')[0]=='N':
                    tr1=tramef_N[station_to_analyse]
                if pair.split(':')[1]=='Z':
                    tr2=tramef_Z[station_to_analyse]
                elif pair.split(':')[1]=='E':
                    tr2=tramef_E[station_to_analyse]
                elif pair.split(':')[1]=='N':
                    tr2=tramef_N[station_to_analyse]

                trames=np.vstack((tr1,tr2))

                del tr1,tr2    
                daycorr = {}
                ndaycorr = {}
                allcorr = {}
                for filterdb in get_filters(db, all=False):
                    filterid = filterdb.ref
                    daycorr[filterid] = np.zeros(get_maxlag_samples(db,))
                    ndaycorr[filterid] = 0

                for islice, (begin, end) in enumerate(zip(begins, ends)):
                    trame2h = trames[:, begin:end]
                    nfft = next_fast_len(int(trame2h.shape[1]))
                    rmsmat = np.std(trame2h, axis=1)
                    for filterdb in get_filters(db, all=False):
                        filterid = filterdb.ref
                        low = float(filterdb.low)
                        high = float(filterdb.high)
                        rms_threshold = filterdb.rms_threshold

                        # Nfft = int(params.min30)
                        # if params.min30 / 2 % 2 != 0:
                        #     Nfft = params.min30 + 2

#Raph                        trames2hWb = np.zeros((2, int(nfft)), dtype=np.complex)
                        print('Shapes of trame2h ', len(trame2h[0]), len(trame2h[1]), len(trame2h))
                        trames2ht = np.zeros((2, int(len(trame2h[0]))), dtype=np.complex)
                        skip = False
                        for i, component in enumerate(pair.split(':')):#Raph changed to split pair, why is STATION still here?
                            if rmsmat[i] > rms_threshold:
                                cp = cosine_taper(len(trame2h[i]),0.04)
                                trame2h[i] -= trame2h[i].mean()

                                if params.windsorizing == -1:
                                    trame2h[i] = np.sign(trame2h[i])
                                elif params.windsorizing != 0:
                                    indexes = np.where(
                                        np.abs(trame2h[i]) > (params.windsorizing * rmsmat[i]))[0]
                                    # clipping at windsorizing*rms
                                    trame2h[i][indexes] = (trame2h[i][indexes] / np.abs(
                                        trame2h[i][indexes])) * params.windsorizing * rmsmat[i]

#Raph                                trames2hWb[i] = whiten(
#                                    trame2h[i]*cp, nfft, dt, low, high, plot=False)
                                print('Testing shapes', trame2h.shape, len(trame2h[0]))                                
                                trames2ht[i] = nowhiten(trame2h[i]*cp, nfft, dt, low, high, plot=False)
                            else:
#Raph                                trames2hWb[i] = np.zeros(int(nfft))
                                trames2ht[i] = np.zeros(int(len(trame2h[i])))
                                skip = True
                                logging.debug('Slice RMS is smaller (%e) than rms_threshold (%e)!'
                                              % (rmsmat[i], rms_threshold))
                        if not skip:
                            corr = myPCorr(trames2ht, np.ceil(params.maxlag / dt), plot=False)#Raph included trames2hFT
                            tmptime = time.gmtime(basetime + begin /
                                                  params.goal_sampling_rate)
                            thisdate = time.strftime("%Y-%m-%d", tmptime)
                            thistime = time.strftime("%Y-%m-%d %H:%M:%S",
                                                     tmptime)
                            print('the shape of the corr is ', corr.shape)
                            if params.keep_all or params.keep_days:#Raph
                                ccfid = "%s_%s_%s_%s" % (station,
                                                         filterid, components,
                                                         thisdate)
                                if ccfid not in allcorr:
                                    allcorr[ccfid] = {}
                                allcorr[ccfid][thistime] = corr

                            if params.keep_days:
                                if not np.any(np.isnan(corr)) and \
                                        not np.any(np.isinf(corr)):
                                    daycorr[filterid] += corr
                                    ndaycorr[filterid] += 1

                            del corr, thistime, trames2ht

                if params.keep_all:
                    for ccfid in allcorr.keys():
                        export_allcorr(db, ccfid, allcorr[ccfid])

                if params.keep_days:
                    for ccfid in allcorr.keys():
                        station, filterid, components, date = ccfid.split('_')#Raph

                        corrs = np.asarray(list(allcorr[ccfid].values()))
                        corr = stack(db, corrs)

                        thisdate = time.strftime(
                                    "%Y-%m-%d", time.gmtime(basetime))
                        thistime = time.strftime(
                                    "%H_%M", time.gmtime(basetime))
                        add_corr(
                                db, station.replace('.', '_'), station.replace('.', '_'), int(filterid),#Raph
                                thisdate, thistime,  params.min30 /
                                params.goal_sampling_rate,
                                components, corr,
                                params.goal_sampling_rate, day=True,
                                ncorr=corrs.shape[0])
                del trames, daycorr, ndaycorr
            logging.debug("Updating Job")
            update_job(db, goal_day, orig_pair, 'AC', 'D')

            logging.info("Finished processing this pair. It took %.2f seconds" % (time.time() - tt))
        logging.info("Job Finished. It took %.2f seconds" % (time.time() - jt))
    logging.info('*** Finished: Compute AC ***')

if __name__ == "__main__":
    main()    

