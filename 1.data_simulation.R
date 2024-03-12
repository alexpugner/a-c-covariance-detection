#clear variables
rm(list=ls(all=TRUE))
#load in libraries
library(MASS)
library(geepack) 
library(sys)
#start a time counter for measurement
start_time=Sys.time()
# ------------------------------------------------------------------------------
# ---------------------------------- settings-----------------------------------
# ------------------------------------------------------------------------------
#--set random seed
set.seed(2101)
# set the number of replication per parameter setting
nrep=500 

# set GEE regression error covariance structure 

cmethod='independence'  

#set the sample size

nmz=2000  
ndz=2000  

#set parameters for Additive genetic effects, shared environmental effects and and non-shared environmental effects simulating data
# These have to add up to 1!!
par_as=p_as=sqrt(c(.34))
par_cs=p_cs=sqrt(c(.14)) 
par_es=p_es=sqrt(c(.53)) 

# set sources of A-C covariance
#g in the model: parents genes have an effect on twins behaviour
par_gs=p_gs=sqrt(c(0,0.025,0.05)) 
#twin 1 genes has an effect on twin 2 behvaiour
par_bs=p_bs=sqrt(c(0,0.025,0.05))
# twin one behaviour has an effect on twin 2 behaviour, and vica-versa
par_xs=p_xs=sqrt(c(0,0.025,0.05)) 
#number of potential maximum sets, for the actual analysis it will be less due to constraints
nset=length(par_as)*length(par_cs)*length(par_es)*length(p_gs)*length(p_bs)*length(p_xs)
nset
#
#setting up 2 matrices, one for the set and one for the results
setkeep=matrix(NA,nset,38)    
setkeep
reskeep=matrix(NA,nset*nrep,38)   


# ------ setting up the values based on quantitative genetics principles
ng=50  # number of diallelic loci   
ngp=10  # number of loci comprising polygenic score pgs (npg shoudl be >= 0 and npg should be <= ng).
p_pgs=ngp/ng  # percentage of genetic variance explained by Polygenic Risk Scores
#
p_A=1-p_pgs # Precentage of genetic variance not explained by polygenic risk scores
#setting up model values
VA1=p_A; VP=p_pgs;VC=1; VE=1
#
ii=0  # count sets
jj=0  # count from 1:(nset*nrep) 
for (par_a in par_as){
  for(par_c in par_cs){
    for (par_e in par_es){
      for (par_g in par_gs) {
        for (par_b in par_bs) {
          for (par_x in par_xs) {
            if (par_b != 0 && par_x !=0) next
              ii=ii+1 # counter per set
              # fill the data frame with set parameters
              setkeep[ii,1:11]=c(nmz,ndz, par_a, par_c, par_e, par_g, par_b, par_x, p_pgs, p_A,ii)
              #
              for (irep in 1:nrep) {  ###### terminates at line +/- 352
                jj=jj+1 #counter per test
                #
                # start ordinary simulation. We set Major and Minor Allele frequencies to .5.
                #
                pal=.5     # maf all GVs same maf ... Minor Allele Freqs (maf)
                qal=1-pal  # major allele freq  
                varGV=2*pal*qal  # variance of genetic variants per SNP
                bs=rep(1,ngp)  
                NMZ=nmz 
                NDZ=ndz   # number of mz and dz twin pairs
                nfam=nmz+ndz # number of families
                
                #--------------------------------------------------------------- 
                # --------------real data simulation starts---------------------
                #---------------------------------------------------------------
                
                
                # simulate parental alleles assuming random mating and linkage equilibrium
                #---------------------------------------------------------------
                #-------------setting up matrices and arrays for----------------
                #---------------------------------------------------------------
                #mother alleles
                am=array(0,c(nfam,ng,2))
                #father alleles
                af=array(0,c(nfam,ng,2)) 
                # mother genotype
                gm=matrix(0,nfam,ng)
                # father genotype
                gf=matrix(0,nfam,ng) 
                # mother pgs
                Am=matrix(0,nfam,1)
                # father pgs
                Af=matrix(0,nfam,1) 
                # mother polygenic risk scores
                pgsm=matrix(0,nfam,1) 
                # father polygenic risk scores
                pgsf=matrix(0,nfam,1) 
                # twin 1 alleles transmitted from mother and father
                at1=array(0,c(nfam,ng,2))
                # twin 1 alleles transmitted from mother and father
                at2=array(0,c(nfam,ng,2))
                # t1 alleles not transmitted from mother and father
                ant1=array(0,c(nfam,ng,2))  
                #twin 2 alleles not transmitted from mother and father
                ant2=array(0,c(nfam,ng,2))  
                #genotypes twin 1
                g1=matrix(0,nfam,ng)
                # genotypes twin 2
                g2=matrix(0,nfam,ng)
                # pgs for twin 1 based on transmitted and nontransmitted alleles (1 = t, 2 = nt)
                A1=matrix(0,nfam,2)
                # pgs for twin 2 based on transmitted and nontransmitted alleles (1 = t, 2 = nt)
                A2=matrix(0,nfam,2) 
                # twin 1 polygenic scores 1=transmitted 2=not transmitted 
                pgs1=matrix(0,nfam,2)
                # twin 2 polygenic scores 1=transmitted 2=not transmitted
                pgs2=matrix(0,nfam,2) 
                
                #---------------------------------------------------------------
                #---------------simulate parental alleles-----------------------
                #---------------------------------------------------------------
                #
                for (i in 1:ng) {
                  am[,i,1] = sample(c(0,1),nfam,replace=T,prob=c(pal,qal))  # mother allele 1
                  am[,i,2] = sample(c(0,1),nfam,replace=T,prob=c(pal,qal)) # mother allele 2
                  af[,i,1] = sample(c(0,1),nfam,replace=T,prob=c(pal,qal)) # father allele 1
                  af[,i,2] = sample(c(0,1),nfam,replace=T,prob=c(pal,qal)) # father allele 2
                }
                
                #---------------------------------------------------------------
                #---------offspring inheritance based on parents alleles--------
                #---------------------------------------------------------------
                
                for (i in 1:ng) {
                  
                  # offspring (twins) inherits 2 alleles per parent
                  # alleles encoded as 1 and 2
                  # maternal allele
                  mt1=sample(c(1,2),nfam,replace=T,prob=c(.5,.5)) 
                  # paternal alleles
                  ft1=sample(c(1,2),nfam,replace=T,prob=c(.5,.5))
                  #calculate non-inherited (non-transmitted) alleles
                  mnt1=3-1*mt1  # 1->2 2-> 1  
                  fnt1=3-1*ft1 # 1->2 2-> 1  
                  
                  # same procedure for twin 2
                  
                  mt2=sample(c(1,2),nfam,replace=T,prob=c(.5,.5)) 
                  ft2=sample(c(1,2),nfam,replace=T,prob=c(.5,.5)) 
                  mnt2=3-1*mt2  
                  fnt2=3-1*ft2
                  # populating matrices for offspring
                  for (k in 1:nfam) {
                    at1[k,i,1] = am[k,i,mt1[k]]   # transm
                    at1[k,i,2] = af[k,i,ft1[k]]   # transm
                    ant1[k,i,1]= am[k,i,mnt1[k]]  # not transm
                    ant1[k,i,2]= af[k,i,fnt1[k]]  # not transm
                    at2[k,i,1] = am[k,i,mt2[k]]   # transm
                    at2[k,i,2] = af[k,i,ft2[k]]   # transm
                    ant2[k,i,1]= am[k,i,mnt2[k]]  # not transm
                    ant2[k,i,2]= af[k,i,fnt2[k]]  # not transm
                  } # nfam
                } # ng
                #
                # A scores representing total additive genetic effects
                #
                for (i in 1:ng) {
                  Am[,1]=Am[,1]+(am[,i,1]+am[,i,2])
                  Af[,1]=Af[,1]+(af[,i,1]+af[,i,2])
                  A1[,1]=A1[,1]+(at1[,i,1]+at1[,i,2])
                  A1[,2]=A1[,2]+(ant1[,i,1]+ant1[,i,2])
                  A2[,1]=A2[,1]+(at2[,i,1]+at2[,i,2])
                  A2[,2]=A2[,2]+(ant2[,i,1]+ant2[,i,2])
                }
                #
                # Calculating polygenic scores based on first ngp genetic variants
                # weighted by bs
                #
                for (i in 1:ngp) {
                  pgsm[,1]=pgsm[,1]+bs[i]*(am[,i,1]+am[,i,2])
                  pgsf[,1]=pgsf[,1]+bs[i]*(af[,i,1]+af[,i,2])
                  pgs1[,1]=pgs1[,1]+bs[i]*(at1[,i,1]+at1[,i,2])
                  pgs1[,2]=pgs1[,2]+(ant1[,i,1]+ant1[,i,2])
                  pgs2[,1]=pgs2[,1]+bs[i]*(at2[,i,1]+at2[,i,2])
                  pgs2[,2]=pgs2[,2]+(ant2[,i,1]+ant2[,i,2])
                }
                #----------------------------------------------------------------------
                #------Standardize variables so the variance adds up to 1-------------
                #---------------------------------------------------------------------
                
                C1_=scale(rnorm(nfam,0,1))  # shared Environment of the twins  C1_ variance = 1 
                #
                E1=scale(rnorm(nfam,0,1))  #   unshared Environment  variance = 1
                E2=scale(rnorm(nfam,0,1))  #   unshared Environment variance = 1
                
                Am[,1]=scale(Am[,1])   # polygenic score mother standardized
                Af[,1]=scale(Af[,1])   # polygenic score father standardized
                A1[,1]=scale(A1[,1])   # polygenic score tw1 standardized (transmitted alleles)
                A1[,2]=scale(A1[,2])   # polygenic score tw1 standardized (not transmitted alleles)
                A2[,1]=scale(A2[,1])   # polygenic score tw2 standardized (transmitted alleles)
                A2[,2]=scale(A2[,2])   # polygenic score tw2 standardized (not transmitted alleles)
                #
                # scale polygenic scores
                #
                pgsm[,1]=scale(pgsm[,1])   # polygenic score mother standardized
                pgsf[,1]=scale(pgsf[,1])   # polygenic score father standardized
                pgs1[,1]=scale(pgs1[,1])   # polygenic score tw1 standardized (transmitted alleles)
                pgs1[,2]=scale(pgs1[,2])   # polygenic score tw1 standardized (not transmitted alleles)
                pgs2[,1]=scale(pgs2[,1])   # polygenic score tw2 standardized (transmitted alleles)
                pgs2[,2]=scale(pgs2[,2])   # polygenic score tw2 standardized (not transmitted alleles)
                #
                # A C E 
                #          DZ
                C11= C1_ # twin 1 
                C12= C1_ # twin 2  ....  
                #           A              C          E     AC cov        AC cov         AC cov (not due to parents)     
                ph1_=par_a*A1[,1] + par_c*C11 + par_e*E1 + par_g*Am[,1] + par_g*Af[,1] + par_b*A2[,1] # offspring 1 pheno 
                ph2_=par_a*A2[,1] + par_c*C12 + par_e*E2 + par_g*Am[,1] + par_g*Af[,1] + par_b*A1[,1] # offspring 1 pheno# offspring 2 pheno 
                # siblinf interaction: mutual influence
                ph1=ph1_+par_x*ph2_
                ph2=ph2_+par_x*ph1_
                #
                # 1:ndz  
                #
                phdatdz=matrix(0,ndz,2)
                phdatdz[,1]=ph1[1:ndz]
                phdatdz[,2]=ph2[1:ndz]
                cov(phdatdz)
                cor(phdatdz)
                #
                # MZ       MZ
                #
                C11= C1_ # twin 1 
                C12= C1_ # twin 2  ....  
                #          A1 and A1                                                           A1 and A1
                ph1_=par_a*A1[,1] + par_c*C11 + par_e*E1 + par_g*Am[,1] + par_g*Af[,1] + par_b*A1[,1] # offspring 1 pheno 
                ph2_=par_a*A1[,1] + par_c*C12 + par_e*E2 + par_g*Am[,1] + par_g*Af[,1] + par_b*A1[,1]  # offspring 1 pheno# offspring 2 pheno 
                # interaction
                ph1=ph1_+par_x*ph2_
                ph2=ph2_+par_x*ph1_
                #
                # ndz+1 : nfam .... second half
                phdatmz=matrix(0,nmz,2)
                phdatmz[,1]=ph1[(ndz+1):nfam]  # ndz+1 to nfam MZs 
                phdatmz[,2]=ph2[(ndz+1):nfam]  # ndz+1 to nfam MZs
                #
                #cov(phdatmz)
                #cor(phdatmz)
                #
                # phdatmz ... add polygenic scores
                #                 mother     father      dz1 t, nt   dz2 t,nt
                pgdatdz=cbind(pgsm[1:ndz,],pgsf[1:ndz,],pgs1[1:ndz,],pgs2[1:ndz,])
                #                 mother            father                  mz1 t, nt =   mz2 t, nt (duplicate 1)       
                pgdatmz=cbind(pgsm[(1+ndz):nfam,],pgsf[(1+ndz):nfam,],pgs1[(1+ndz):nfam,],pgs1[(1+ndz):nfam,])
                #
                # stochastic simulated data 
                phdatmz=as.data.frame(cbind(phdatmz,pgdatmz))
                phdatdz=as.data.frame(cbind(phdatdz,pgdatdz))
                # add sum and mean
                colnames(phdatdz) = colnames(phdatmz) =  vnames1 = 
                  c('pht1','pht2',
                    'pgsm','pgsf','pgst1','pgsnt1','pgst2','pgsnt2')
                adddz=cbind(phdatdz$pgsm+phdatdz$pgsf, (phdatdz$pgst1+phdatdz$pgst2)/2)
                colnames(adddz) = c('pgsmf','mpgst')
                phdatdz = cbind(phdatdz, adddz)
                addmz=cbind(phdatmz$pgsm+phdatmz$pgsf, (phdatmz$pgst1+phdatmz$pgst2)/2)
                colnames(addmz) = c('pgsmf','mpgst')
                phdatmz = cbind(phdatmz, addmz)
                #
                c(1,2,3,4,5,7,9,10) -> i1
                
                #
                apply(phdatmz,2,var) 
                round(cor(phdatmz),3) 
                #
                apply(phdatdz,2,var) 
                round(cor(phdatdz),3) 
                #
                # [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst" 
                # Organize data in long format simulated data
                #
                phdatmzL = matrix(0,nmz*2,8)
                phdatmzL[,2]=c(c(1:nmz),c(1:nmz))
                phdatmzL[,3]=c(phdatmz$pht1,phdatmz$pht2)
                phdatmzL[,4]=c(phdatmz$pgsm,phdatmz$pgsm)
                phdatmzL[,5]=c(phdatmz$pgsf,phdatmz$pgsf)
                phdatmzL[,6]=c(phdatmz$pgst1,phdatmz$pgst2)
                phdatmzL[,7]=c(phdatmz$pgsmf,phdatmz$pgsmf)
                phdatmzL[,8]=c(phdatmz$mpgst,phdatmz$mpgst)
                #
                ix_ = sort.int(phdatmzL[,2], index.return=T)
                phdatmzL = phdatmzL[ix_$ix,]
                colnames(phdatmzL)=c("zyg","famnr","ph","pgsm","pgsf","pgst","pgsmf","mpgst")
                phdatmzL=as.data.frame(phdatmzL)
                #
                # Long format simulated data 
                phdatdzL = matrix(1,ndz*2,8)
                phdatdzL[,2]=nmz + c(c(1:ndz),c(1:ndz))
                phdatdzL[,3]=c(phdatdz$pht1,phdatdz$pht2)
                phdatdzL[,4]=c(phdatdz$pgsm,phdatdz$pgsm)
                phdatdzL[,5]=c(phdatdz$pgsf,phdatdz$pgsf)
                phdatdzL[,6]=c(phdatdz$pgst1,phdatdz$pgst2)
                phdatdzL[,7]=c(phdatdz$pgsmf,phdatdz$pgsmf)
                phdatdzL[,8]=c(phdatdz$mpgst,phdatdz$mpgst)
                #
                ix_ = sort.int(phdatdzL[,2], index.return=T)
                phdatdzL = phdatdzL[ix_$ix,]
                colnames(phdatdzL)=c("zyg","famnr","ph","pgsm","pgsf","pgst","pgsmf","mpgst")
                phdatdzL=as.data.frame(phdatdzL)
                #
                phdatL=rbind(phdatmzL, phdatdzL)
                #
                # simulated stochastically
                #                wide        wide     long      long      long mz+dz
                # data sets are  phdatdz and phdatdz, phdatdzL, phdatmzL, phdatL
                #                wide          wide        long       long        long mz+dz
                # simulated exactly
                # data sets are  phdatdz_e and phdatdz_e, phdatdzL_e, phdatmzL_e, phdatL_e
                #              pheno t1 pheno t2 mother    father   twin 1  nt twin1  twin2   nt twin2  m+f prs  mean twin prs
                # colnames [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
                # 
                # --------------------------------------------------- end data sim
                #
                # regression analyses. based on simulated data (not exact) ...
                #
                # DZ twin 1 only ... just pgsmf
                #
                M0dz=lm(pht1~pgst1, data=phdatdz)#)$coefficients       	#  single twin dz 1    ... pheno on pgs ...just regression
                M1dz=lm(pht1~pgsmf+pgst1, data=phdatdz)#)$coefficients 	#  single twin dz 1    ... pheno on pgs + m&f pgs test of cov(AC)
                #M2dz=lm(pht1~mpgst+pgst1, data=phdatdz)#)$coefficients 	#    single twin dz 1   ... pheno on pgs + twin mean pgs test of cov AC ... not important
                #M3dz=lm(pht1~pgsmf+mpgst+pgst1, data=phdatdz)#)$coefficients #     single twin dz 1 ... pheno on pgs, m&f pgs + twin mean pgs test of cov AC
                #
                # 2 dz twins regression using gee 
                #
                geeM0dzL=geeglm(ph~pgst, id=famnr, corstr=cmethod, data=phdatdzL)#)$coefficients #  on pgs   
                geeM1dzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod, data=phdatdzL)#)$coefficients #   m+f pgs + twin pgs 
                geeM2dzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL)#)$coefficients  #   mean twin pgs + twin pgs
                geeM3dzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL)#)$coefficients  #    # mf+mean+twin
                #
                geeM0mzdzL=geeglm(ph~pgst, id=famnr, corstr=cmethod,data=phdatL)#)$coefficients #          same sequence of analysis
                geeM1mzdzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod,data=phdatL)#)$coefficients #   
                geeM2mzdzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL)#)$coefficients  #   
                geeM3mzdzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL)#)$coefficients  #  
                
                R2M0mzdzL=summary(lm(ph~pgst,data=phdatL))$r.squared #          same sequence of analysis
                R2M1mzdzL=summary(lm(ph~pgsmf+pgst,data=phdatL))$r.squared #  
                R2M2mzdzL=summary(lm(ph~mpgst+pgst,data=phdatL))$r.squared  #  
                R2M3mzdzL=summary(lm(ph~pgsmf+mpgst+pgst,data=phdatL))$r.squared  #     # mf+mean+twin
                #
                R2M0dzL=summary(lm(ph~pgst,data=phdatdzL))$r.squared #          same sequence of analysis
                R2M1dzL=summary(lm(ph~pgsmf+pgst,data=phdatdzL))$r.squared #  
                R2M2dzL=summary(lm(ph~mpgst+pgst,data=phdatdzL))$r.squared  #  
                R2M3dzL=summary(lm(ph~pgsmf+mpgst+pgst,data=phdatdzL))$r.squared  #  
                #
                # reskeep=matrix(NA,nrep,30)   # to keep results
                # get power exact power pgsmf, mpgst, pgsmf+mpgst: 3 tests ...  dz pairs,   dz + mz pairs
                # get power exact power pgsmf: 1 tests ...  dz singles
                #
                # dz 1 pgsmf test ... test of mpgst does not apply given 1 dz
                # model M1dz=lm(pht1~pgsmf+pgst1, data=phdatdz) ... 1dz test of pgsmf
                tmp=summary(M1dz)$coefficients
                M1dzest=tmp[2,1]; M1dzse=tmp[2,2]; M1dzp=tmp[2,4]
                reskeep[jj,1:3]=c(M1dzest,M1dzse,M1dzp) # est st and power
                #
                # dz 1+2 test power
                # model geeM1dzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod, data=phdatdzL) #   
                # pgsmf
                tmp=summary(geeM1dzL)$coefficients
                geeM1dzest=tmp[2,1]; geeM1dzse=tmp[2,2]; geeM1dzp=tmp[2,4]
                reskeep[jj,4:6]=c(geeM1dzest,geeM1dzse,geeM1dzp) # est st and power
                # dz 1+2 test power
                #geeM2dzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL)#)$coefficients  #   
                # mpgst
                tmp=summary(geeM2dzL)$coefficients
                geeM2dzest=tmp[2,1]; geeM2dzse=tmp[2,2]; geeM2dzp=tmp[2,4]
                reskeep[jj,7:9]=c(geeM2dzest,geeM2dzse,geeM2dzp) # est st and power
                #
                # geeM3dzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL)#)$coefficients  #   
                # pgsmf in presence of mfpgs
                tmp=summary(geeM3dzL)$coefficients
                geeM3dzest1=tmp[2,1]; geeM3dzse1=tmp[2,2]; geeM3dzp1=tmp[2,4]; 
                reskeep[jj,10:12]=c(geeM3dzest1,geeM3dzse1,geeM3dzp1) # est st and power
                # mpgst
                geeM3dzest2=tmp[3,1]; geeM3dzse2=tmp[3,2]; geeM3dzp2=tmp[3,4]
                reskeep[jj,13:15]=c(geeM3dzest2,geeM3dzse2,geeM3dzp2) # est st and power
                #
                # mz dz 1+2 test power
                # geeM1mzdzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients #
                # pgsmf
                tmp=summary(geeM1mzdzL)$coefficients
                geeM1mzdzest=tmp[2,1]; geeM1mzdzse=tmp[2,2]; geeM1mzdzp=tmp[2,4];
                reskeep[jj,16:18]=c(geeM1mzdzest,geeM1mzdzse,geeM1mzdzp) # est st and power
                # dz 1+2 test power
                # egeeM2mzdzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients  #   
                # mpgst
                tmp=summary(geeM2mzdzL)$coefficients
                geeM2mzdzest=tmp[2,1]; geeM2mzdzse=tmp[2,2]; geeM2mzdzp=tmp[2,4]
                reskeep[jj,19:21]=c(geeM2mzdzest,geeM2mzdzse,geeM2mzdzp) # est st and power
                #
                # geeM3mzdzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients  # #   
                # pgsmf in presence of mfpgs
                tmp=summary(geeM3mzdzL)$coefficients
                geeM3mzdzest1=tmp[2,1]; geeM3mzdzse1=tmp[2,2]; geeM3mzdzp1=tmp[2,4]; 
                reskeep[jj,22:24]=c(geeM3mzdzest1,geeM3mzdzse1,geeM3mzdzp1) # est st and power
                # mpgst
                geeM3mzdzest2=tmp[3,1]; geeM3mzdzse2=tmp[3,2]; geeM3mzdzp2=tmp[3,4];
                reskeep[jj,25:27]=c(geeM3mzdzest2,geeM3mzdzse2,geeM3mzdzp2) # est st and power
                #
                reskeep[jj,28:35]=c(R2M0mzdzL,R2M1mzdzL,R2M2mzdzL,R2M3mzdzL,R2M0dzL,R2M1dzL,R2M2dzL,R2M3dzL)
                reskeep[jj,36:38]=c(ii,irep,jj)
              } # irep         line 352 about
            }}}}}} # the factorial design
setkeep
colnames(reskeep) = c('dz1mfpgs_e','dz1mfpgs_s','dz1mfpgs1_p',
                      'dz2mfpgs_e','dz2mfpgs_s','dz2mfpgs_p',
                      'dz2mpgsm_e','dz2mpgsm_s','dz2mpgsm_p',
                      'dz2mfpgs_m_e','dz2mfpgs_m_s','dz2mfpgs_m_p',
                      'dz2mpgst_mf_e','dz2mpgst_mf_s','dz2mpgst_mf_p',
                      'mzdzmfpgs_e','mzdzmfpgs_s','mzdzmfpgs_p',
                      'mzdzmpgsm_e','mzdzmpgsm_s','mzdzmpgsm_p',
                      'mzdzmfpgs_m_e','mzdzmfpgs_m_s','mzdzmfpgs_m_p',
                      'mzdzmpgst_mf_e','mzdzmpgst_mf_s','mzdzmpgst_mf_p',
                      'R2M0mzdzL','R2M1mzdzL','R2M2mzdzL','R2M3mzdzL','R2M0dzL',
                      'R2M1dzL','R2M2dzL','R2M3dzL','set','irep','tcounter')
reskeep
write_csv(as.data.frame(setkeep), '/Users/macszerez.com/Desktop/portfolio/conor/data/setkeep.csv')
write_csv(as.data.frame(reskeep), '/Users/macszerez.com/Desktop/portfolio/conor/data/reskeep.csv')
#------------------deleting NA rows (caused by the conditional statement in the for loop) 
#------------------and creating 3 different matrices with standard errors, estimates and p-values)
end_time=Sys.time()
time_elapsed=end_time-start_time
time_elapsed




