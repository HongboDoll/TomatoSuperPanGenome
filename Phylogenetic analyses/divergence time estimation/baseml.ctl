
      seqfile = single_copy_gene_family_pep_4_tree.phylip * sequence data file name
      outfile = subsititution_rate        * main result file
     treefile = tomato_single_copy_phyml_tree_rooted.nwk_baseml  * tree structure file name

        noisy = 3   * 0,1,2,3: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

        model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV
		    * 8:UNREST, 9:REVu; 10:UNRESTu
        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

    fix_kappa = 0
        kappa = 2   * initial or given kappa

    fix_alpha = 0 
        alpha = 0.5  * initial or given alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates

      fix_rho = 1  
          rho = 0  * initial or given rho,   0:no correlation
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 
		   
        clock = 1   * 0: no clock, unrooted tree, 1: clock, rooted tree
        nhomo = 1   * 0 & 1: homogeneous, 2: kappa's, 3: N1, 4: N2
        getSE = 1   * 0: don't want them, 1: want S.E.s of estimates
	method = 0 * Optimization method 0: simultaneous;1: one branch a time
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
