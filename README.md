# _FECGSYN_ toolbox, version 1.2, March 2017

Open-source platform for reproducible NI-FECG research

<center>
<img src="http://fernandoandreotti.github.io/fecgsyn/pages/images/French_fecgsyn.png" width="200">
</center>

* Documentation available at: [www.fecgsyn.com](http://www.fecgsyn.com)
* Large dataset of simulated signals (FECGSYNDB) [available at Physionet](http://physionet.org/physiobank/database/fecgsyndb/)
* See latest release under [https://github.com/fernandoandreotti/fecgsyn/releases/](https://github.com/fernandoandreotti/fecgsyn/releases/)


## Authors

_FECGSYN_ is the product of a collaboration between the Department of Engineering Science, University of Oxford (DES-OX), the Institute of Biomedical Engineering, TU Dresden (IBMT-TUD) and the Biomedical Engineering Faculty at the Technion Israel Institute of Technology (BME-IIT). The authors are:
- Joachim Behar (DES-OX,BME-IIT)
- Fernando Andreotti (DES-OX,IBMT-TUD)
- Julien Oster (DES-OX)
- Sebastian Zaunseder (IBMT-TUD) 
- Gari D. Clifford (DES-OX)


Contacts: joachim.behar@oxfordalumni.org & fernando.andreotti@eng.ox.ac.uk


## History

FECGSYN is built upon the work from McSharry et al. [1] and Sameni et al. [2] 
The original code from McSharry et al. is available in MATLAB and in 
C on PhysioNet. The code developed by Sameni et al. 
is part of the OSET toolbox, also available online in MATLAB.
Links to these work are available at: 
http://physionet.incor.usp.br/physiotools/ipmcode/

1. McSharry, Patrick E and Clifford, Gari D and Tarassenko, Lionel and Smith, Leonard A.
A dynamical model for generating synthetic electrocardiogram signals. IEEE Transactions
on Biomedical Engineering,  50(3) 2003.

2. Sameni, Reza, et al. Multichannel ECG and noise modeling: application to
maternal and foetal ECG signals. EURASIP Journal on Advances in Signal Processing
2007 (2007).

## References


When using the _FECGSYN_ please refer to:

Behar, Joachim, Fernando Andreotti, Sebastian Zaunseder, Qiao Li, Julien Oster, and Gari D Clifford. 2014. 
"An ECG Model for Simulating Maternal-Foetal Activity Mixtures on Abdominal ECG Recordings." _Physiol Meas_ **35(8)**, pp.1537-50, 2014.

Andreotti F., Behar J., Zaunseder S.,Oster J. and Clifford G D., An Open-Source Framework for Stress-Testing Non-Invasive Foetal ECG Extraction Algorithms. _Physiol Meas_ **37(5)**, pp. 627-648, 2016.



## License


Released under the GNU General Public License

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

## Disclaimer

This toolbox makes use of several other pre-existing open source algorithms listed below:

- **ECGSYN: A realistic ECG waveform generator**, by Dr. Patrick McSharry and Gari D. Clifford,  [https://www.physionet.org/physiotools/ecgsyn/](available here) (licensed under GNU GPL 2.0)
- **Open Source ECG Toolbox (OSET)**, v1.0, by Dr. Reza Sameni, [http://oset.ir/](available here) (licensed under GNU GPL 2.0)
- **FastICA for Matlab**, v2.5, by  Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen [http://research.ics.aalto.fi/ica/fastica/](available here) (licensed under GNU GPL 2.0)
- **JadeR**, by Jean-Francois Cardoso, [http://perso.telecom-paristech.fr/~cardoso](available here) (BSD license)
- **FECG-ESN toolbox**, v1.0, Dr. Joachim Behar, [http://joachimbehar.comuv.com](available here)  (licensed under GNU GPL 2.0)
- **ESN learning toolbox**, v1.0, by H. Jaeger (Fraunhofer IAIS), [http://reservoir-computing.org/software](available here) (unlicensed)
- **QRS Detection with Pan-Tompkins algorithm**, by Daniel Wedekind, [https://github.com/danielwedekind/qrsdetector](available here)  (licensed under GNU GPL 2.0)
- **arrow.m**, by Dr. Erik A. Johnson, [https://uk.mathworks.com/matlabcentral/fileexchange/278-arrow](available here), (BSD license)
- **fwhm**, v1.2, by Patrick Egan, [http://uk.mathworks.com/matlabcentral/fileexchange/10590-fwhm](available here), (BSD license)
- **pcorr2**, by Peter Rydesäter, [https://uk.mathworks.com/matlabcentral/fileexchange/4012-prcorr2-10-times-faster-correlation-coef](available here) (BSD license)

Not provided with package, ocasionally required, see [http://fernandoandreotti.github.io/fecgsyn/pages/install.html](install instructions):
- **WFDB Toolbox for MATLAB and Octave**, v.0.9.9, by Dr. Ikaro Silva, [https://www.physionet.org/physiotools/matlab/wfdb-app-matlab/](available here) (licensed under GNU GPL 2.0)
