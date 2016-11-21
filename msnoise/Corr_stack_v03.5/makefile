# call : make "name of program"
FFLAG1=-mcmodel=medium
FFLAG1=
FFLAG2=
libsac=/usr/local/sac/lib/sacio.a
MYBIN=.
MYSRC=.

pcc5g: $(MYSRC)/pcc5g.f
	gfortran $(FFLAG1) $(FFLAG2) $(MYSRC)/pcc5g.f -o $(MYBIN)/pcc5g $(libsac)

tf_pws0: $(MYSRC)/tf_pws0.f
	gfortran $(FFLAG1) $(FFLAG2) $(MYSRC)/tf_pws0.f -o $(MYBIN)/tf_pws0 $(libsac)

norm_tavg_white: $(MYSRC)/norm_tavg_white.f
	gfortran $(FFLAG1) $(FFLAG2) $(MYSRC)/norm_tavg_white.f -o $(MYBIN)/norm_tavg_white $(libsac)
