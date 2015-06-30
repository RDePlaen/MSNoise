""" This code is responsible for the computation of the cross-correlation
functions.

This script will group *jobs* marked "T"odo in the database by day and process
them using the following scheme. As soon as one day is selected, the
corresponding jobs are marked "I"n Progress in the database. This allows
running several instances of this script in parallel.


Waveform Preprocessing
~~~~~~~~~~~~~~~~~~~~~~
Pairs are first split and a station list is created. The database is then
queried to get file paths. For each station, all files potentially containing
data for the day are opened. The traces are then merged and splitted, to obtain
the most continuous chunks possible. The different chunks are then demeaned,
tapered and merged again to a 1-day long trace. If shorter than 1-day, the
trace is padded with zeros. If longer, it is cut to match the start/end of the
day.

Each 1-day long trace is then low-passed (at ``preprocess_lowpass`` Hz),
high-passed (at ``preprocess_highpass`` Hz), then decimated/downsampled.
Decimation/Downsampling are configurable (``resampling_method``) and users are
advised testing both. One advantage of Downsampling over Decimation is that
it is able to downsample the data by any factor, not only integer factors.

.. warning::
    For an unknown reason, the PAZ-correction has disappeard from the
    current sqlvolution on GitHub: CHECK!


Processing
~~~~~~~~~~

Once all traces are preprocessed, station pairs are processed sequentially.
If a component different from *ZZ* is to be computed, the traces are first
rotated. This supposes the user has provided the station coordinates in the
*station* table. The rotation is computed for Radial and Transverse components:

.. code-block:: python

    R = tramef_N * np.cos(cplAz * np.pi / 180.) + tramef_E * np.sin(cplAz * np.pi / 180.)
    T = tramef_N * np.sin(cplAz * np.pi / 180.) - tramef_E * np.cos(cplAz * np.pi / 180.)

Then, for each ``corr_duration`` window in the signal, and for each filter
configured in the database, the traces are clipped to ``windsorizing`` times
the RMS and then whitened (see :ref:`whiten`) between the frequency bounds.
When both traces are ready, the cross-correlation function is computed
(see :ref:`mycorr`). The function returned contains data for time lags 
corresponding to ``maxlag`` in the acausal (negative lags) and causal
(positive lags) parts.

If configured (setting ``keep_all`` to 'Y'), each ``corr_duration`` CCF is
saved to the hard disk. By default, the ``keep_days`` setting is set to True
and so "N = 1 day / corr_duration" CCF are stacked and saved to the hard disk
in the STACKS/001_DAYS folder.

Once done, each job is marked "D"one in the database.

To run this script:

.. code-block:: sh

    $ msnoise compute_cc
"""

import time
import calendar
import sys
from obspy import read
from obspy.core import utcdatetime
from scikits.samplerate import resample
import scipy.fftpack

from api import *
from myCorr import myCorr
#from whiten import whiten


def preprocess(db, stations, comps, goal_day, params, tramef_Z, tramef_E=np.array([]), tramef_N=np.array([])):
    datafilesZ = {}
    datafilesE = {}
    datafilesN = {}

    for station in stations:
        datafilesZ[station] = []
        datafilesE[station] = []
        datafilesN[station] = []
        net, sta = station.split('.')
        gd = datetime.datetime.strptime(goal_day, '%Y-%m-%d')
        files = get_data_availability(
            db, net=net, sta=sta, starttime=gd, endtime=gd)
        for file in files:
            comp = file.comp
            fullpath = os.path.join(file.path, file.file)
            if comp[-1] == 'Z':
                datafilesZ[station].append(fullpath)
            elif comp[-1] == 'E':
                datafilesE[station].append(fullpath)
            elif comp[-1] == 'N':
                datafilesN[station].append(fullpath)

    j = 0
    for istation, station in enumerate(stations):
        for comp in comps:
            files = eval("datafiles%s['%s']" % (comp, station))
            if len(files) != 0:
                logging.debug("%s.%s Reading %i Files" %
                              (station, comp, len(files)))
                stream = Stream()
                for file in sorted(files):
                    st = read(file, dytpe=np.float)
                    for tr in st:
                        tr.data = tr.data.astype(np.float)
                    stream += st
                    del st

                logging.debug("Checking sample alignment")
                for i, trace in enumerate(stream):
                    stream[i] = check_and_phase_shift(trace)

                stream.sort()
                logging.debug("Checking Gaps")
                if len(getGaps(stream)) > 0:
                    max_gap = 10
                    only_too_long = False
                    while getGaps(stream) and not only_too_long:
                        too_long = 0
                        gaps = getGaps(stream)
                        for gap in gaps:
                            if int(gap[-1]) <= max_gap:
                                stream[gap[0]] = stream[gap[0]].__add__(stream[gap[1]], method=0,
                                                                        fill_value="interpolate")
                                stream.remove(stream[gap[1]])
                                break
                            else:
                                too_long += 1
                        if too_long == len(gaps):
                            only_too_long = True

                taper_length = 20.0  #seconds
                for trace in stream:
                    if trace.stats.npts < 4 * taper_length * trace.stats.sampling_rate:
                        trace.data = np.zeros(trace.stats.npts)
                    else:
                        trace.detrend(type="demean")
                        trace.detrend(type="linear")
                        taper_1s = taper_length * float(trace.stats.sampling_rate) / trace.stats.npts
                        cp = cosTaper(trace.stats.npts, taper_1s)
                        trace.data *= cp
                stream.merge(method=0, fill_value=0.0)

                logging.debug("%s.%s Slicing Stream to %s:%s" % (station, comp, utcdatetime.UTCDateTime(
                    goal_day.replace('-', '')), utcdatetime.UTCDateTime(
                    goal_day.replace('-', '')) + params.goal_duration - stream[0].stats.delta))
                stream[0].trim(utcdatetime.UTCDateTime(goal_day.replace('-', '')), utcdatetime.UTCDateTime(
                    goal_day.replace('-', '')) + params.goal_duration - stream[0].stats.delta, pad=True, fill_value=0.0,
                               nearest_sample=False)
                trace = stream[0]

                logging.debug(
                    "%s.%s Highpass at %.2f Hz" % (station, comp, params.preprocess_highpass))
                trace.filter("highpass", freq=params.preprocess_highpass, zerophase=True)

                if trace.stats.sampling_rate != params.goal_sampling_rate:
                    logging.debug(
                        "%s.%s Lowpass at %.2f Hz" % (station, comp, params.preprocess_lowpass))
                    trace.filter("lowpass", freq=params.preprocess_lowpass, zerophase=True)

                    if params.resampling_method == "Resample":
                        logging.debug("%s.%s Downsample to %.1f Hz" %
                                      (station, comp, params.goal_sampling_rate))
                        trace.data = resample(
                            trace.data, params.goal_sampling_rate / trace.stats.sampling_rate, 'sinc_fastest')

                    elif params.resampling_method == "Decimate":
                        logging.debug("%s.%s Decimate by a factor of %i" %
                                      (station, comp, params.decimation_factor))
                        trace.data = trace.data[::params.decimation_factor]
                    trace.stats.sampling_rate = params.goal_sampling_rate

                year, month, day, hourf, minf, secf, wday, yday, isdst = trace.stats.starttime.utctimetuple()

                if j == 0:
                    t = time.strptime("%04i:%02i:%02i:%02i:%02i:%02i" %
                                      (year, month, day, hourf, minf, secf), "%Y:%m:%d:%H:%M:%S")
                    basetime = calendar.timegm(t)

                if len(trace.data) % 2 != 0:
                    trace.data = np.append(trace.data, 0.)

                if comp == "Z":
                    tramef_Z[istation] = trace.data
                elif comp == "E":
                    tramef_E[istation] = trace.data
                elif comp == "N":
                    tramef_N[istation] = trace.data

                del trace, stream
    if len(tramef_E) != 0:
        return basetime, tramef_Z, tramef_E, tramef_N
    else:
        return basetime, tramef_Z


class Params():
    pass


def main():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('*** Starting: Compute AC ***')

    # Connection to the DB
    db = connect()
    #rule out absence of filters
    if len(get_filters(db, all=False)) == 0:
        print "NO FILTER!!"
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
    params.components_to_compute = ['Z', 'E', 'N']
    ####Select here only the comps in get_data_availability
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
    #for station_unique in stations_to_analyse:
    while is_next_job(db, jobtype='AC'):
        jobs = get_next_job(db, jobtype='AC')
        #logging.info("Working on station %s" %(station_unique))

        stations = []
        refs = []

        #go through job to make job array with stations
        for job in jobs:
            refs.append(job.ref)
            #pairs.append(job.pair)#find a way to pair comps   /!\
            netsta = job.pair  #just 1 station in the cell?
            stations.append(netsta)
            #stations.append(netstacomp2)
            goal_day = job.day

        stations = np.unique(stations)#only 1 station in it
        #print "only this station(s) here: %s, ref= %s" %(stations, refs)
        logging.info("New AC Job: %s (%i pairs with %i stations)" %
                     (goal_day, len(pairs)*len(stations), len(stations)))
        jt = time.time()

        xlen = int(params.goal_duration * params.goal_sampling_rate)


        comps = ['Z', 'E', 'N']
        tramef_Z = np.zeros((len(stations), xlen))
        tramef_E = np.zeros((len(stations), xlen))
        tramef_N = np.zeros((len(stations), xlen))
        basetime, tramef_Z, tramef_E, tramef_N = preprocess(db, stations, comps, goal_day, params, tramef_Z, tramef_E, tramef_N)# preprocessing
        #print type(tramef_E)

        # comps = ['Z']
        # tramef_Z = np.zeros((len(stations), xlen))
        # basetime, tramef_Z = preprocess(db, stations, comps, goal_day, params, tramef_Z)
        # print type(tramef_Z)

        dt = 1. / params.goal_sampling_rate
        # Calculate the number of slices

        slices = int(params.goal_duration * params.goal_sampling_rate / params.min30)
        begins = []
        ends = []
        i = 0
        while i <= (params.goal_duration - params.min30 / params.goal_sampling_rate):
            begins.append(int(i * params.goal_sampling_rate))
            ends.append(int(i * params.goal_sampling_rate + params.min30))
            i += int(params.min30 / params.goal_sampling_rate * (1.0 - params.overlap))
        slices = len(begins)


        # ##########################################################################################################

        for station in stations:
            orig_pair = station
            for pair in pairs:
                print "Processing pair %s for station %s" %(pair, station)
                logging.info('Processing pair: %s' % pair)
                #print type(tramef_Z)
                #print tramef_Z.keys()
                tt = time.time()
                comp1,comp2=pair.split(':')
                components=comp1+comp2
                ### load trames
                #assign trames according to pair

                station_to_analyse=np.where(stations==station)
                #print "you are looking for ",station_to_analyse



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
                #print "tr1 est de type ",type(tr1)

                trames=np.vstack((tr1,tr2))

#                print tr1
#                print tr2

                del tr1,tr2

                ## islice
                daycorr = {}
                ndaycorr = {}
                allcorr = {}
                for filterdb in get_filters(db, all=False):
                    filterid = filterdb.ref
                    daycorr[filterid] = np.zeros(get_maxlag_samples(db,))
                    ndaycorr[filterid] = 0

                for islice, (begin, end) in enumerate(zip(begins, ends)):
                    #print "Progress: %#2d/%2d"% (islice+1,slices)
                    trame2h = trames[:, begin:end]

                    rmsmat = np.std(np.abs(trame2h), axis=1)
                    for filterdb in get_filters(db, all=False):
                        filterid = filterdb.ref
                        low = float(filterdb.low)
                        high = float(filterdb.high)
                        rms_threshold = filterdb.rms_threshold

                        Nfft = int(params.min30)
                        if params.min30 / 2 % 2 != 0:
                            Nfft = params.min30 + 2

                        #trames2hWb = np.zeros((2, int(Nfft)), dtype=np.complex)
                        trames2hFT = np.zeros((2, int(Nfft)), dtype=np.complex)
                        skip = False
                        for i, comp in enumerate(pair.split(':')):#### fix this loop, work on each comp of the pair
                            #print i, comp, station, goal_day
                            if rmsmat[i] > rms_threshold:
                                cp = cosTaper(len(trame2h[i]),0.04)
                                trame2h[i] -= trame2h[i].mean()

                                if params.windsorizing == -1:
                                    trame2h[i] = np.sign(trame2h[i])
                                elif params.windsorizing != 0:
                                    indexes = np.where(
                                        np.abs(trame2h[i]) > (params.windsorizing * rmsmat[i]))[0]
                                    # clipping at windsorizing*rms
                                    trame2h[i][indexes] = (trame2h[i][indexes] / np.abs(
                                        trame2h[i][indexes])) * params.windsorizing * rmsmat[i]

                                #trames2hWb[i] = whiten(trame2h[i]*cp, Nfft, dt, low, high, plot=False)
                                trames2hFT[i]=scipy.fftpack.fft(trame2h[i]*cp, Nfft)
                            else:
                                #trames2hWb[i] = np.zeros(int(Nfft))
                                trames2hFT[i] = np.zeros(int(Nfft))
                                skip = True
                                logging.debug('Slice is Zeros!')

                        if not skip:
                            corr = myCorr(trames2hFT, np.ceil(params.maxlag / dt), plot=False)
                            tmptime = time.gmtime(basetime + begin /
                                                  params.goal_sampling_rate)
                            thisdate = time.strftime("%Y-%m-%d", tmptime)
                            thistime = time.strftime("%Y-%m-%d %H:%M:%S",
                                                     tmptime)
                            if params.keep_all:
                                ccfid = "%s_%s_%s_%s_%s" % (station,
                                                         filterid, components,
                                                         thisdate)
                                if ccfid not in allcorr:
                                    allcorr[ccfid] = {}
                                allcorr[ccfid][thistime] = corr

                            if params.keep_days:
                                #print "KEEP DAYS!"
                                if not np.any(np.isnan(corr)) and \
                                        not np.any(np.isinf(corr)):
                                    daycorr[filterid] += corr
                                    ndaycorr[filterid] += 1

                            del corr, thistime, trames2hFT
                        else:
                            print "NOOOOOOOOOOO! Zeros!"

                if params.keep_all:
                    for ccfid in allcorr.keys():
                        export_allcorr(db, ccfid, allcorr[ccfid])

                if params.keep_days:
                    try:
                        for filterdb in get_filters(db, all=False):
                            filterid = filterdb.ref
                            corr = daycorr[filterid]
                            ncorr = ndaycorr[filterid]
                            if ncorr > 0:
                                logging.debug(
                                    "Saving daily CCF for filter %02i, comp %s (stack of %02i CCF)" % (filterid, components, ncorr))

                                thisdate = time.strftime(
                                    "%Y-%m-%d", time.gmtime(basetime))
                                thistime = time.strftime(
                                    "%H_%M", time.gmtime(basetime))
                                stationpair="%s_%s" %(station,station)
                                add_corr(
                                    db, station.replace('.', '_'),#################
                                    station.replace('.', '_'), filterid,
                                    thisdate, thistime,  params.min30 /
                                    params.goal_sampling_rate,
                                    components, corr,
                                    params.goal_sampling_rate, day=True,
                                    ncorr=ncorr)
                                del corr, ncorr
                    except Exception as e:
                            logging.debug(str(e))
                del trames, daycorr, ndaycorr
            logging.debug("Updating Job")
            update_job(db, goal_day, orig_pair, 'AC', 'D')

            logging.info("Finished processing this station. It took %.2f seconds" %
                              (time.time() - tt))
        logging.info("Job Finished. It took %.2f seconds" % (time.time() - jt))
logging.info('*** Finished: Compute AC ***')

if __name__ == "__main__":
    main()
