NI-FECG simulator toolbox, version 1.0, February 2014
====
Released under the GNU General Public License


Copyright (C) 2014 
Intelligent Patient Monitoring Group - Oxford 2014
Contacts: joachim.behar@eng.ox.ac.uk & fernando.andreotti@mailbox.tu-dresden.de

Description
FECGSYN is a realistic non-invasive foetal ECG (NI-FECG) generator that 
uses the Gaussian ECG model originally introduced by McSharry et al. 
The toolbox generates synthetic NI-FECG mixtures considering various 
user-defined settings, e.g. noise sources, heart rate and heart rate 
variability, rotation of the maternal and foetal heart axes due to 
respiration, foetus movement, contractions, ectopic beats and multiple 
pregnancy. Any number of electrodes can be freely placed on the maternal 
abdomen. The synthetic ECG simulator is a good tool for modelling realistic 
FECG-MECG mixtures and specific events such as abrupt heart rate increase, 
in order to benchmark signal processing algorithms on realistic data and 
for scenarios that resemble important clinical events.

FECGSYN is fruit of the collaboration between the Department of Engineering 
Science, University of Oxford (DES-OX) and the Institute of Biomedical Engineering, 
TU Dresden (IBMT-TUD). The authors are Joachim Behar (DES-OX), Fernando Andreotti 
(IBMT-TUD), Julien Oster (DES-OX), Sebastian Zaunseder (IBMT-TUD) and 
Gari Clifford (DES-OX). 

History
FECGSYN is built upon the work from McSharry et al. and Sameni et al. 
The original code from McSharry et al. is available in MATLAB and in 
C on PhysioNet. The code developed by Sameni et al. 
is part of the OSET toolbox, also available online in MATLAB.
Links to these work are available at: 
http://physionet.incor.usp.br/physiotools/ipmcode/



This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

For further reference:

Behar, Joachim, Fernando Andreotti, Sebastian Zaunseder, Qiao Li, Julien Oster, and Gari D Clifford. 2014. 
?An ECG Model for Simulating Maternal-Foetal Activity Mixtures on Abdominal ECG Recordings.? 
Physiol Meas (Accepted for publication - Focus issue: Noninvasive Fetal ECG).

Bibtex

@article{Behar_Andreotti_Zaunseder_Li_Oster_Clifford_2014, 
title={An ECG Model for Simulating Maternal-Foetal Activity Mixtures on Abdominal ECG Recordings}, 
author={Behar, Joachim and Andreotti, Fernando and Zaunseder, Sebastian and Li, Qiao and Oster, Julien and Clifford, Gari D}, 
number={Focus issue: Noninvasive Fetal ECG}, 
journal={Physiol Meas (Accepted for publication)}, 
year={2014}}
