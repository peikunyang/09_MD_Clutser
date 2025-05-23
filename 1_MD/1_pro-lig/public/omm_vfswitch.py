from openmm.unit import *
from openmm import *
from openmm.app import *

def vfswitch(system, psf, inputs):
    r_on = inputs.r_on
    r_off = inputs.r_off

    # custom nonbonded force for force-switch
    chknbfix = False
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nonbonded = force
        if isinstance(force, CustomNonbondedForce) and force.getNumTabulatedFunctions() == 2:
            nbfix     = force
            chknbfix  = True

    # vfswitch
    vfswitch = CustomNonbondedForce('step(Ron-r)*(ccnba*tr6*tr6-ccnbb*tr6+ccnbb*onoff3-ccnba*onoff6) \
                                     +step(r-Ron)*step(Roff-r)*(cr12*rjunk6 - cr6*rjunk3) \
                                     -step(r-Ron)*step(Ron-r)*(cr12*rjunk6 - cr6*rjunk3); \
                                     cr6  = ccnbb*ofdif3*rjunk3; \
                                     cr12 = ccnba*ofdif6*rjunk6; \
                                     rjunk3 = r3-recof3; \
                                     rjunk6 = tr6-recof6; \
                                     r3 = r1*tr2; \
                                     r1 = sqrt(tr2); \
                                     tr6 = tr2 * tr2 * tr2; \
                                     tr2 = 1.0/s2; \
                                     s2 = r*r; \
                                     ccnbb = 4.0*epsilon*sigma^6; \
                                     ccnba = 4.0*epsilon*sigma^12; \
                                     sigma = sigma1+sigma2; \
                                     epsilon = epsilon1*epsilon2; \
                                     onoff3 = recof3/on3; \
                                     onoff6 = recof6/on6; \
                                     ofdif3 = off3/(off3 - on3); \
                                     ofdif6 = off6/(off6 - on6); \
                                     recof3 = 1.0/off3; \
                                     on6 = on3*on3; \
                                     on3 = c2onnb*Ron; \
                                     recof6 = 1.0/off6; \
                                     off6 = off3*off3; \
                                     off3 = c2ofnb*Roff; \
                                     c2ofnb = Roff*Roff; \
                                     c2onnb = Ron*Ron; \
                                     Ron  = %f; \
                                     Roff = %f;' % (r_on, r_off) )
    vfswitch.addPerParticleParameter('sigma')
    vfswitch.addPerParticleParameter('epsilon')
    vfswitch.setNonbondedMethod(vfswitch.CutoffPeriodic)
    vfswitch.setCutoffDistance(nonbonded.getCutoffDistance())
    for i in range(nonbonded.getNumParticles()):
        chg, sig, eps = nonbonded.getParticleParameters(i)
        nonbonded.setParticleParameters(i, chg, 0.0, 0.0) # zero-out LJ
        sig = sig*0.5
        eps = eps**0.5
        vfswitch.addParticle([sig, eps])
    for i in range(nonbonded.getNumExceptions()):
        atom1, atom2 = nonbonded.getExceptionParameters(i)[:2]
        vfswitch.addExclusion(atom1, atom2)
    vfswitch.setForceGroup(psf.NONBONDED_FORCE_GROUP)
    system.addForce(vfswitch)

    # vfswitch14
    vfswitch14 = CustomBondForce('step(Ron-r)*(ccnba*tr6*tr6-ccnbb*tr6+ccnbb*onoff3-ccnba*onoff6) \
                                  +step(r-Ron)*step(Roff-r)*(cr12*rjunk6 - cr6*rjunk3) \
                                  -step(r-Ron)*step(Ron-r)*(cr12*rjunk6 - cr6*rjunk3); \
                                  cr6  = ccnbb*ofdif3*rjunk3; \
                                  cr12 = ccnba*ofdif6*rjunk6; \
                                  rjunk3 = r3-recof3; \
                                  rjunk6 = tr6-recof6; \
                                  r3 = r1*tr2; \
                                  r1 = sqrt(tr2); \
                                  tr6 = tr2 * tr2 * tr2; \
                                  tr2 = 1.0/s2; \
                                  s2 = r*r; \
                                  ccnbb = 4.0*epsilon*sigma^6; \
                                  ccnba = 4.0*epsilon*sigma^12; \
                                  onoff3 = recof3/on3; \
                                  onoff6 = recof6/on6; \
                                  ofdif3 = off3/(off3 - on3); \
                                  ofdif6 = off6/(off6 - on6); \
                                  recof3 = 1.0/off3; \
                                  on6 = on3*on3; \
                                  on3 = c2onnb*Ron; \
                                  recof6 = 1.0/off6; \
                                  off6 = off3*off3; \
                                  off3 = c2ofnb*Roff; \
                                  c2ofnb = Roff*Roff; \
                                  c2onnb = Ron*Ron; \
                                  Ron  = %f; \
                                  Roff = %f;' % (r_on, r_off) )
    vfswitch14.addPerBondParameter('sigma')
    vfswitch14.addPerBondParameter('epsilon')
    for i in range(nonbonded.getNumExceptions()):
        atom1, atom2, chg, sig, eps = nonbonded.getExceptionParameters(i)
        nonbonded.setExceptionParameters(i, atom1, atom2, chg, 0.0, 0.0) # zero-out LJ14
        vfswitch14.addBond(atom1, atom2, [sig, eps])
    system.addForce(vfswitch14)

    # vfswitch_NBFIX
    if chknbfix:
        nbfix.setEnergyFunction('step(Ron-r)*(ccnba*tr6*tr6-ccnbb*tr6+ccnbb*onoff3-ccnba*onoff6) \
                                 +step(r-Ron)*step(Roff-r)*(cr12*rjunk6 - cr6*rjunk3) \
                                 -step(r-Ron)*step(Ron-r)*(cr12*rjunk6 - cr6*rjunk3); \
                                 cr6  = ccnbb*ofdif3*rjunk3; \
                                 cr12 = ccnba*ofdif6*rjunk6; \
                                 rjunk3 = r3-recof3; \
                                 rjunk6 = tr6-recof6; \
                                 r3 = r1*tr2; \
                                 r1 = sqrt(tr2); \
                                 tr6 = tr2 * tr2 * tr2; \
                                 tr2 = 1.0/s2; \
                                 s2 = r*r; \
                                 ccnbb = bcoef(type1, type2); \
                                 ccnba = acoef(type1, type2)^2; \
                                 onoff3 = recof3/on3; \
                                 onoff6 = recof6/on6; \
                                 ofdif3 = off3/(off3 - on3); \
                                 ofdif6 = off6/(off6 - on6); \
                                 recof3 = 1.0/off3; \
                                 on6 = on3*on3; \
                                 on3 = c2onnb*Ron; \
                                 recof6 = 1.0/off6; \
                                 off6 = off3*off3; \
                                 off3 = c2ofnb*Roff; \
                                 c2ofnb = Roff*Roff; \
                                 c2onnb = Ron*Ron; \
                                 Ron  = %f; \
                                 Roff = %f;' % (r_on, r_off) )

        # turn off long range correction (OpenMM Issues: #2353)
        nbfix.setUseLongRangeCorrection(False)

    return system
