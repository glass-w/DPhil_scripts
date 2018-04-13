##########################
## A macro to select the three chains
## in the beta3 trimer model built
###########################

#### Chains ####

set b3_trimer_model.CHNA {(resid 1 to 166)}
atomselect macro CHNA ${b3_trimer_model.CHNA}

set b3_trimer_model.CHNB {(resid 167 to 332)}
atomselect macro CHNB ${b3_trimer_model.CHNB}

set b3_trimer_model.CHNC {(resid 333 to 498)}
atomselect macro CHNC ${b3_trimer_model.CHNC}
