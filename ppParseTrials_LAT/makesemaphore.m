function [ns6directory_sem] = makesemaphore(ns6directory)
%Gives the appropriate semaphore for each arrayID
%Make sure the arrayID and sessions to look at match 

ns6directory_sem=ns6directory;

%Beaker
ns6directory_sem(1).Semaphore=[100 0 100 0 100 0 1 0 1 0 1 0];
ns6directory_sem(3).Semaphore=[100 0 100 0 100 0 1 0 1 0 1 0];
ns6directory_sem(8).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(12).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(14).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(18).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(20).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(22).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(26).Semaphore=[100 100 100 1 0 1 0 1 0];

%Scooter
ns6directory_sem(28).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(30).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(32).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(34).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(36).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(38).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(40).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(42).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(44).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(46).Semaphore=[0 1 0 1 0];


end

