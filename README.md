# functional_programming_demo
A brief demo to illustrate the benefits of fp.

This is part of a presentation about functional programming and its perks. The demo provides a naive Matlab implementation of anomaly detection algorithm HOTSAX (https://ieeexplore.ieee.org/document/1565683/), and a timeseries from QT database (https://www.physionet.org/physiobank/database/qtdb/) onto which anomaly detection is performed.

apply_func_windows is the central function which takes care of looping through the timeseries, windowing it and applying a function to each window. Once that function is done and dusted (written once, used happily ever after), we can concentrate on the actual problem at hand and write an implementation for HOTSAX with a few lines of re-usable code. 
