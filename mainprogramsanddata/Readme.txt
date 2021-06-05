To replicate Tables 2 and 7 and the ID diagnostic tests, one needs the final data set, wdset.dat.
Put the programs and the data sets in the same directory.

The rest of the programs are simulation programs that do not require data.

Table 2:             sumstat.g, sumstata.g, pcb.g wdset.dat (The first two programs calculate basic summary statistics, and the third calculates autocorrelation coefficients.)
Table 3:             gammaK.g        (You will need to do two runs, one with dofe set to 0 and another with dofe set to 1. This flag is near the top of the program. You will also need to set the sample size = 1500.)
Table 4:             gammaA.g        (You will need to do two runs, one with dofe set to 0 and another with dofe set to 1. This flag is near the top of the program. You will also need to set the sample size = 1500.)
Table 5:             gammaK.g, gammaA.g  (Here you need to set dofe=1, but set dofe=0 to get the results for the IV estimators. You will need to set the sample size to 4000.)

                     For Tables 3-5, you are going to have to Monte Carlo a bootstrap. This takes a very long time. The program EWbootMC.g is designed to do exactly one Monte Carlo trial and then spit out ones and zeros according
                     to whether the bootstrapped test rejected or not. You will get 6 numbers:

                     test for the EW estimate of beta
                     test for the EW estimate of alpha1
                     test for the IV estimate of beta
                     test for the IV estimate of alpha1
                     test for the AB estimate of beta
                     test for the AB estimate of alpha1

                     You need the input file, params.txt, which contains five numbers that specify the following:

                     1 for fixed effects, 0 for levels
                     1 for the third moment estimstor, 2 for the fourth moment estimator, 3 for the fifth moment estimator, 4 for the sixth moment estimator, etc.
                     1 for a samples size of 1500, 2 for 4000
                     1 for the capital stock DGP, 2 for the assets DGP
                     The trial number, which goes from 1 to 1000.

                     We did this on the Amazon Elastic Compute Cloud as follows. You will need to run 22 individual Monte Carlos of 1000 trials each, but because you can do many trials at once on different machines, it goes relatively
                     quickly. Each Monte Carlo corresponds to a column of Tables 3-5. For each one, compile EWbootMC.g. Write a script that
                     makes the file params.txt with all of the right inputs and then runs one Monte Carlo trial with the Gauss run-time module. Then make the script update the trial number and repeat, preferably not all on the same
                     machine. You can do several hundred of these at once, in which you set a different (sequential) trial number for each one. Collect all of the output files, read them in, and average them.

                     Note: to replicate the results exactly, you will have to run one Monte Carlo trial at a time, because the program resets the seed with each run of the program. You will get close if you just run it all at once,
                     (set the variable mcount=1000), but be prepared to wait for about 7 weeks for some of these to finish on a top of the line Xeon chip, circa late 2011.

Table 6:             Everything, including instructions, is in the directory filesfortable6
                    (See also the original files sent to us by Almeida, Campello, and Galvao in ACGunmodified.)


id tests:            idtest1.g, cum34.g, wdset.dat

                     You will need to run each of these programs four times (fe=0,asset=0), (fe=1,asset=0), (fe=0,asset=1), (fe=1,asset=1). You'll see where to change these at the top of the program.


Table 7:             search.g, ewboot.g, wdset.dat

                     You will have to run search.g four times (fe=0,asset=0), (fe=1,asset=0), (fe=0,asset=1), (fe=1,asset=1)
                     Note that for ewboot.g, you are going to have to run one estimator at a time, for each of the above four cases, making a total of 12 runs.  The program prints out bootstrapped p-values
                     in square brackets underneath the coefficient estimates.

                     You will also need to compute the Arellano Bond and IV estimates. The dataset abdset.dat has first-differenced variables for the AB estimation. The programs are ab.g and simpleiv.g.

Figure 1:            gammaKskew.g

Figure 2:            gammaKhet(run twice, dofe=0 and dofe=1), gammaKfe, gammaKcorre.g

Figure 3:            midtest.g, mcum34.g

                     You will need to run each of these programs four times (dofe=0,ncount=1500), (dofe=1,ncount=1500), (dofe=0,ncount=4000), (dofe=1,ncount=4000). You'll see where to change these at the top of the program.
