          seed = -1
       seqfile = single_copy_gene_family_pep_4_tree.phylip
      treefile = tomato_single_copy_phyml_tree_rooted.nwk_mcmctree
       outfile = out.BV

         ndata = 1
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 3    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
*       RootAge = '<1.0'  * safe constraint on root age, used if no fossil for root.

         model = 7    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 7:REV(GTR)
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1  * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 1 2.26 * gammaDir prior for rate for genes (1 / subsititution_rate)
  sigma2_gamma = 1 4.5  * gammaDir prior for sigma^2     (for clock=2 or 3)

      finetune = 1: 0.06  0.5  0.006  0.12 0.4 * times, musigma2, rates, mixing, paras, FossilErr

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 500000
      sampfreq = 100
       nsample = 20000

*** Note: Make your window wider (100 columns) before running the program.
