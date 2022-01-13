% STARTUP  some project- and Alex-specific MATLAB settings for SBT.
% Run this before doing anything.

% compulsory
addpath('utils')

% some of alex's generic startup stuff
com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',false)
format long g
format compact
set(0,'showHiddenHandles','on');        % eg lines linesmoothing
% restore pre-R2014b look:
set(groot,'DefaultFigureGraphicsSmoothing','off')   % kill slow antialiasing
set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                                   .75 0 .75; .75 .75 0; .25 .25 .25])
set(groot,'DefaultFigureColormap',jet(256))
set(groot,'defaultLegendAutoUpdate','off')

maxNumCompThreads(16);   % speeds up some multithreaded matlab ops, on 8-core.
