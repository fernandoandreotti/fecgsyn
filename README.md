[![license](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](./LICENSE)

# _FECGSYN_ toolbox, version 1.3-alpha, August 2019

### Open-source platform for reproducible NI-FECG research

<p align="center"><img src="http://fernandoandreotti.github.io/fecgsyn/pages/images/French_fecgsyn.png" width="200"></p>

* Documentation available at: [www.fecgsyn.com](http://www.fecgsyn.com)
* Large dataset of simulated signals (FECGSYNDB) [available at Physionet](http://physionet.org/physiobank/database/fecgsyndb/)
* See latest release under [https://github.com/fernandoandreotti/fecgsyn/releases/](https://github.com/fernandoandreotti/fecgsyn/releases/)


## Authors

_FECGSYN_ is the product of a collaboration between the Department of Engineering Science, University of Oxford (DES-OX), the Institute of Biomedical Engineering, TU Dresden (IBMT-TUD), the Department of Electrical and Electronic Engineering, University of Melbourne (EEE-UOM) and the Biomedical Engineering Faculty at the Technion Israel Institute of Technology (BME-IIT). The authors are:
- Joachim Behar
- Fernando Andreotti
- Julien Oster
- Sebastian Zaunseder
- Gari D. Clifford
- Emerson Keenan
- Chandan Karmakar
- Marimuthu Palaniswami


Contacts: jbehar@technion.ac.il & fernando.andreotti@eng.ox.ac.uk


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


When using _FECGSYN_ please reference at least one of the following articles:

1. Behar, Joachim, Fernando Andreotti, Sebastian Zaunseder, Qiao Li, Julien Oster, and Gari D Clifford. 2014. "An ECG Model for Simulating Maternal-Foetal Activity Mixtures on Abdominal ECG Recordings." _Physiol Meas_ **35(8)**, pp.1537-50, 2014.

and/or

2. Andreotti F., Behar J., Zaunseder S.,Oster J. and Clifford G D., An Open-Source Framework for Stress-Testing Non-Invasive Foetal ECG Extraction Algorithms. _Physiol Meas_ **37(5)**, pp. 627-648, 2016.

If you are using _FECGSYN's_ asymmetric volume conductor modeling capability, please reference the following article:

3. Keenan E., Karmakar C K. and Palaniswami M., The effects of asymmetric volume conductor modeling on non-invasive fetal ECG extraction. _Physiol Meas_ **39(10)**, pp. 105013, 2018.

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

- _ECGSYN: A realistic ECG waveform generator_, by Dr. Patrick McSharry and Gari D. Clifford,  [available here](https://www.physionet.org/physiotools/ecgsyn/) (licensed under GNU GPL 2.0)
- _Open Source ECG Toolbox (OSET)_, v1.0, by Dr. Reza Sameni, [available here](http://oset.ir/) (licensed under GNU GPL 2.0)
- _FastICA for Matlab_, v2.5, by  Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen [available here](http://research.ics.aalto.fi/ica/fastica/) (licensed under GNU GPL 2.0)
- _JadeR_, by Jean-Francois Cardoso, [available here](http://perso.telecom-paristech.fr/~cardoso) (BSD license)
- _FECG-ESN toolbox_, v1.0, Dr. Joachim Behar, [available here](http://joachimbehar.comuv.com)  (licensed under GNU GPL 2.0)
- _ESN learning toolbox_, v1.0, by H. Jaeger (Fraunhofer IAIS), [available here](http://reservoir-computing.org/software) (unlicensed)
- _QRS Detection with Pan-Tompkins algorithm_, by Daniel Wedekind, [available here](https://github.com/danielwedekind/qrsdetector)  (licensed under GNU GPL 2.0)
- _arrow.m_, by Dr. Erik A. Johnson, [available here](https://uk.mathworks.com/matlabcentral/fileexchange/278-arrow), (BSD license)
- _fwhm.m_, v1.2, by Patrick Egan, [available here](http://uk.mathworks.com/matlabcentral/fileexchange/10590-fwhm), (BSD license)
- _pcorr2.m_, by Peter Rydesäter, [available here](https://uk.mathworks.com/matlabcentral/fileexchange/4012-prcorr2-10-times-faster-correlation-coef) (BSD license)
- _FieldTrip: The MATLAB toolbox for MEG and EEG analysis_, v20190828, by Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen, [available here](https://github.com/fieldtrip/fieldtrip) (licensed under GNU GPL 2.0)
- _Iso2Mesh_, v1.9.0-1, Qianqian Fang, [available here](https://github.com/fangq/iso2mesh) (licensed under GNU GPL 2.0)

Not provided with package, ocasionally required, see [install instructions](http://fernandoandreotti.github.io/fecgsyn/pages/install.html):
- **WFDB Toolbox for MATLAB and Octave**, v.0.9.9, by Dr. Ikaro Silva, [available here](https://www.physionet.org/physiotools/matlab/wfdb-app-matlab/) (licensed under GNU GPL 2.0)
- **Pre-processed anatomic models**, by Emerson Keenan, [available here](https://github.com/emersonkeenan/fecgsyn-anatomic-models) (licensed under CC BY-NC-SA 2.0 FR)
