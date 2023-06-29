function [ns6directory_sem] = makesemaphore(ns6directory)
%Gives the appropriate semaphore for each arrayID

ns6directory_sem=ns6directory;

%Beaker
ns6directory_sem(4).Semaphore=[100 0 100 0 100 0 1 0 1 0 1 0];
ns6directory_sem(9).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(13).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(15).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(19).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(21).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(23).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(25).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(27).Semaphore=[100 100 100 1 0 1 0 1 0];

%Scooter
ns6directory_sem(29).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(31).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(33).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(35).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(37).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(39).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(41).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(43).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(45).Semaphore=[100 100 100 1 0 1 0 1 0];
ns6directory_sem(47).Semaphore=[0 1 0 1 0];


end

