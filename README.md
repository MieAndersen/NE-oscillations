# NE-oscillations
Code to reproduce figures in "Memory enhancing properties of sleep depend on the oscillatory amplitude of norepinephrine" by Kjaerby, Andersen, Hauglund, Untiet, Dall, Sigurdsson, Ding, Feng, Li, Weikop, Hirase, and Nedergaard (2022).

The code is ordered into folders with "scripts" containing scripts used for analysis and "functions" containing custom functions that are used in the scripts.

Additionally, to load fiber photometry data we use a custom function provided by Tucker Davis Technologies called "TDTbin2mat" that you can find in their [repository](https://github.com/tdtneuro/TDTMatlabSDK). To load EEG and EMG data obtained using Sleepscore (ViewPoint, France) we use a custom function provided by ViewPoint Behavior Technology called "loadEXP" - this function van be obtain upon request ([www.viewpoint.fr](www.viewpoint.fr)).
