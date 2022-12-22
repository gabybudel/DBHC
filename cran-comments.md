## Resubmission of archived package
This is a resubmission of an archived package. In this version I have:

* Removed the assumption that class of matrix has length 1 (the error for which the package was archived)
* Made smaller example execute on test
* Added test directory with simple test for errors


## Test environments
* local Windows 10 install, R 4.2.2
* Windows Server 2022, R-devel, 64 bit (R-hub)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (R-hub)
* Fedora Linux, R-devel, clang, gfortran (R-hub)

## R CMD check results
* There were no ERRORs, WARNINGs, or NOTEs on local Windows install and Windows 
Server 2022.
* There were no ERRORs or WARNINGs on Ubuntu and Fedora Linux, but 2 notes:
* Possibly misspelled words in DESCRIPTION:
  DBHC (10:6, 10:27)
  HMMs (3:49, 9:6, 11:55) 
  Answer: spelling is correct.
* checking examples ... [5s/19s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
            user system elapsed
  hmm.clust 1.363  0.028   5.457
On Fedora Linux exceeded max runtime of 5s by a few ms. The example really is 
small, so should be fine.
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Answer: manual is fine.
* Note that it is a resubmission of an archived package on CRAN. Said errors are resolved now.
