import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import colossus.cosmology.cosmology as cosmo
from colossus.halo import  mass_so
import seaborn as sns
from scipy.stats import kde
import matplotlib as mpl
from scipy.signal import savgol_filter as svg
cosmol=cosmo.setCosmology('planck15')    #cosmology for colossus. By default hmf uses planck15
cosmo.setCurrent(cosmol)



plt.rcParams['font.size']=45
plt.rcParams['figure.figsize']=(16,16)
#mpl.rcParams['xtick.minor.visible']=True
plt.rcParams['axes.linewidth']= 3.
plt.rcParams['axes.titlepad'] = 20
plt.rcParams['axes.linewidth']=5
plt.rcParams['xtick.major.size'] =15
plt.rcParams['ytick.major.size'] =15
plt.rcParams['xtick.minor.size'] =10
plt.rcParams['ytick.minor.size'] =10
plt.rcParams['xtick.major.width'] =5
plt.rcParams['ytick.major.width'] =5
plt.rcParams['xtick.minor.width'] =5
plt.rcParams['ytick.minor.width'] =5
plt.rcParams['axes.titlepad'] = 20 

def make_savgol(y,yup,ydown,window):
    
    y = svg(y,window,1)
    yup = svg(yup,window,1)
    ydown = svg(ydown,window,1)

    return y,yup,ydown
    

def make_catalog(df,redshift):
    
    centrals_index = df.groupby('GrpID')['SubDM'].transform(max) == df['SubDM']

    centrals = df[centrals_index]
    sat = df[np.logical_not(centrals_index)]
    strings = ['ParentDM','SubDM','SatSM','Re']
#for s in strings:
#    centrals.loc[:,s]=  centrals.loc[:,s]/0.7

    centrals['Rhalo'] = mass_so.M_to_R(centrals.SubDM.apply(lambda x: 10**x), z=redshift,mdef='200c')

#maskSFR = centrals['SFR'].values==0
#cenrals['SFR'][maskSFR] = 

    for s in strings:
        centrals[s]=  centrals[s] -np.log10(cosmol.h)
        df[s] = df[s] -np.log10(cosmol.h)
        sat[s]= sat[s] -np.log10(cosmol.h)
    
    centrals['Rhalo'] = np.log10(centrals['Rhalo'] -np.log10(cosmol.h))
    centrals['Re'] = np.log10(centrals['Re'] - np.log10(cosmol.h))

    return centrals


def make_percentiles(rhalo,re):
    
    rhalobins = np.arange(1.5,3,0.2)
    med = []
    up =[]
    low=[]
    
    for i in range(len(rhalobins)-1):
        
        try:
            mask = np.ma.masked_inside(rhalo, rhalobins[i], rhalobins[i+1]).mask
        
            re_ = re[mask]
        
            res = np.percentile(re_,[16,50,84])
            low.append(res[0])
            med.append(res[1])
            up.append(res[2])
        except:
            low.append(np.nan)
            med.append(np.nan)
            up.append(np.nan)
            
        
    return 10**np.array(low), 10**np.array(med),10**np.array(up)
    

if __name__ == '__main__':
'''    
    df0 = pd.read_csv('./IllSnaps/z_0.0.csv')
    df1 = pd.read_csv('./IllSnaps/z_0.9973.csv')
    df2 = pd.read_csv('./IllSnaps/z_2.002.csv')

    cen0 = make_catalog(df0,0.0)
    cen1 = make_catalog(df1,0.9973)
    cen2 = make_catalog(df2,2.002)

    cen0['SSFR'] = cen0['SFR']/(cen0['SatSM'].apply(lambda x: 10**x))*1.e9
    cen1['SSFR'] = cen1['SFR']/(cen1['SatSM'].apply(lambda x: 10**x))*1.e9
    cen2['SSFR'] = cen2['SFR']/(cen2['SatSM'].apply(lambda x: 10**x))*1.e9

    cen0['SSFR'][cen0.SSFR==0] = 1.e-5
    cen1['SSFR'][cen1.SSFR==0] = 1.e-5
    cen2['SSFR'][cen2.SSFR==0] = 1.e-5
    
    cen0 = cen0[(cen0.SatSM>9) ]#& (cen0.SSFR>0)]
    cen1 = cen1[(cen1.SatSM>9) ]#& (cen1.SSFR>0)]
    cen2 = cen2[(cen2.SatSM>9) ]#& (cen2.SSFR>0)]
    
    rhalobins = np.arange(1.5,3,0.2)
    
    ridges = np.array([-2,-2,-2]) #np.array([-0.94,-0.35,0.05]) - 1.
    centrals_cut = cen0[(cen0.SatSM.values>9)]# & (cen0.SSFR>0) ]
    starforming0 = centrals_cut[np.log10(centrals_cut.SSFR.values)>ridges[0] ]
    quiescent0 = centrals_cut[np.log10(centrals_cut.SSFR.values)< ridges[0] ]
    
    low0_SF,med0_SF,up0_SF = make_percentiles(starforming0.Rhalo.values, starforming0.Re.values)
    low0_Q,med0_Q,up0_Q = make_percentiles(quiescent0.Rhalo.values, quiescent0.Re.values)
    
    
    centrals_cut = cen1[(cen1.SatSM.values>9)]# & (cen1.SSFR>0) ]
    starforming1 = centrals_cut[np.log10(centrals_cut.SSFR.values)>ridges[1] ]
    quiescent1 = centrals_cut[np.log10(centrals_cut.SSFR.values)<ridges[1] ]
    
    low1_SF,med1_SF,up1_SF = make_percentiles(starforming1.Rhalo.values, starforming1.Re.values)
    low1_Q,med1_Q,up1_Q = make_percentiles(quiescent1.Rhalo.values, quiescent1.Re.values)
    
    centrals_cut = cen2[(cen2.SatSM.values>9)]# & (cen2.SSFR>0) ]
    starforming2 = centrals_cut[np.log10(centrals_cut.SSFR.values)>ridges[2] ]
    quiescent2 = centrals_cut[np.log10(centrals_cut.SSFR.values)<ridges[2] ]
    
    low2_SF,med2_SF,up2_SF = make_percentiles(starforming2.Rhalo.values, starforming2.Re.values)
    low2_Q,med2_Q,up2_Q = make_percentiles(quiescent2.Rhalo.values, quiescent2.Re.values)
    
    
    fig,(ax1,ax2,ax3) = plt.subplots(3,3,figsize = (32,32),sharey=True,sharex=True)
    
    
    rhalobins = 10**(rhalobins[1:]-0.1)
    ax1[0].plot(rhalobins,med0_SF,color='dodgerblue',lw=5,label='LTGs')
    ax1[0].fill_between(rhalobins, low0_SF, up0_SF, color= 'dodgerblue', alpha=0.4)
    
    upscat = 10**(np.log10(med0_SF) + 0.2)
    lowscat = 10**(np.log10(med0_SF) - 0.2)
    
    ax1[0].plot(rhalobins, upscat, ls=':', color='dodgerblue', alpha=0.5,lw=4)
    ax1[0].plot(rhalobins, lowscat, ls=':', color='dodgerblue', alpha=0.5,lw=4)
    
    #####################################
    
    ax1[1].plot(rhalobins,med1_SF,color='dodgerblue',lw=5)
    ax1[1].fill_between(rhalobins, low1_SF, up1_SF, color= 'dodgerblue', alpha=0.4)

    upscat = 10**(np.log10(med1_SF) + 0.2)
    lowscat = 10**(np.log10(med1_SF) - 0.2)
    
    ax1[1].plot(rhalobins, upscat, ls=':', color='dodgerblue', alpha=0.5,lw=4)
    ax1[1].plot(rhalobins, lowscat, ls=':', color='dodgerblue', alpha=0.5,lw=4)
    
    #########################################
    
    ax1[2].plot(rhalobins,med2_SF,color='dodgerblue',lw=5)
    ax1[2].fill_between(rhalobins, low2_SF, up2_SF, color= 'dodgerblue', alpha=0.4)
   
    upscat = 10**(np.log10(med2_SF) + 0.2)
    lowscat = 10**(np.log10(med2_SF) - 0.2)
    
    ax1[2].plot(rhalobins, upscat, ls=':', color='dodgerblue', alpha=0.5,lw=4)
    ax1[2].plot(rhalobins, lowscat, ls=':', color='dodgerblue', alpha=0.5,lw=4)            
       
      ##################################
    
    
    ax1[0].plot(rhalobins,med0_Q,color='orangered',lw=5,label='ETGs')
    ax1[0].fill_between(rhalobins, low0_Q, up0_Q, color= 'orangered', alpha=0.4)

    upscat = 10**(np.log10(med0_Q) + 0.2)
    lowscat = 10**(np.log10(med0_Q) - 0.2)
    
    ax1[0].plot(rhalobins, upscat, ls=':', color='orangered', alpha=0.5, lw=4)
    ax1[0].plot(rhalobins, lowscat, ls=':', color='orangered', alpha=0.5, lw=4)
    
     #####################################
        
    ax1[1].plot(rhalobins,med1_Q,color='orangered',lw=5)
    ax1[1].fill_between(rhalobins, low1_Q, up1_Q, color= 'orangered', alpha=0.4)

    upscat = 10**(np.log10(med1_Q) + 0.2)
    lowscat = 10**(np.log10(med1_Q) - 0.2)
    
    ax1[1].plot(rhalobins, upscat, ls=':', color='orangered', alpha=0.5, lw=4)
    ax1[1].plot(rhalobins, lowscat, ls=':', color='orangered', alpha=0.5, lw=4)
    
    ###########################
    
    ax1[2].plot(rhalobins,med2_Q,color='orangered',lw=5)
    ax1[2].fill_between(rhalobins, low2_Q, up2_Q, color= 'orangered', alpha=0.4)
    
    upscat = 10**(np.log10(med2_Q) + 0.2)
    lowscat = 10**(np.log10(med2_Q) - 0.2)
    
    ax1[2].plot(rhalobins, upscat, ls=':', color='orangered', alpha=0.5, lw=4)
    ax1[2].plot(rhalobins, lowscat, ls=':', color='orangered', alpha=0.5, lw=4)  
    
    
    ############################################### SAM ########################################



    #################  ssfr selected ##############
    
    z=1
    x0_SF,med0_SF,up0_SF,down0_SF = np.loadtxt('./SAM/ssfr_selected/ltgs/z'+str(z)+'/median.dat', unpack=True)
    x0_Q,med0_Q,up0_Q,down0_Q = np.loadtxt('./SAM/ssfr_selected/etgs/z'+str(z)+'/median.dat', unpack=True)

    m = np.ma.masked_greater(x0_Q,50).mask

    x0_Q = x0_Q[m]
    
    #med0_Q,up0_Q,down0_Q = make_savgol(med0_Q[m],up0_Q[m],down0_Q[m], window=3)

    z=3
    x1_SF,med1_SF,up1_SF,down1_SF = np.loadtxt('./SAM/ssfr_selected/ltgs/z'+str(z)+'/median.dat', unpack=True)
    x1_Q,med1_Q,up1_Q,down1_Q = np.loadtxt('./SAM/ssfr_selected/etgs/z'+str(z)+'/median.dat', unpack=True)

    
    m = np.ma.masked_greater(x1_Q,50).mask
    x1_Q = x1_Q[m]
   # med1_Q,up1_Q,down1_Q = make_savgol(med1_Q[m],up1_Q[m],down1_Q[m],window=5)
    
    
    z=5
    x2_SF,med2_SF,up2_SF,down2_SF = np.loadtxt('./SAM/ssfr_selected/ltgs/z'+str(z)+'/median.dat', unpack=True)
   # x2_Q,med2_Q,up2_Q,down2_Q = np.loadtxt('./SAM/ssfr_selected/etgs/z'+str(z)+'/median.dat', unpack=True)

    m = np.ma.masked_greater(x2_Q,50).mask
    x2_Q = x2_Q[m]
    med2_Q = med2_Q[m]
    up2_Q = up2_Q[m]
    down2_Q = down2_Q[m]

   # med2_Q,up2_Q,down2_Q = make_savgol(med2_Q,up2_Q,down2_Q)

   
    ax2[0].plot(x0_SF,med0_SF, ls='-',  color='dodgerblue',lw=5)
    ax2[0].fill_between(x0_SF,up0_SF,down0_SF, color='dodgerblue',alpha=0.4)
    ax2[0].plot(x0_Q,med0_Q, ls='-',  color='orangered',lw=5)
    ax2[0].fill_between(x0_Q,up0_Q,down0_Q, color='orangered',alpha=0.4)
    
    upSF_ = 10**(np.log10(med0_SF)+0.2)
    lowSF_ = 10**(np.log10(med0_SF)-0.2)
    upQ_ = 10**(np.log10(med0_Q)+0.2)
    lowQ_ = 10**(np.log10(med0_Q)-0.2)
    
    ax2[0].plot(x0_SF,upSF_,color='dodgerblue',ls=':',lw=4)
    ax2[0].plot(x0_SF,lowSF_,color='dodgerblue',ls=':',lw=4)    
    ax2[0].plot(x0_Q,upQ_,color='orangered',ls=':',lw=4)
    ax2[0].plot(x0_Q,lowQ_,color='orangered',ls=':',lw=4)
    
    
    ax2[1].plot(x1_SF,med1_SF, ls='-',  color='dodgerblue',lw=5)
    ax2[1].fill_between(x1_SF,up1_SF,down1_SF, color='dodgerblue',alpha=0.4)
    ax2[1].plot(x1_Q,med1_Q, ls='-',  color='orangered',lw=5)
    ax2[1].fill_between(x1_Q,up1_Q,down1_Q, color='orangered',alpha=0.4)    

    
    upSF_ = 10**(np.log10(med1_SF)+0.2)
    lowSF_ = 10**(np.log10(med1_SF)-0.2)
    upQ_ = 10**(np.log10(med1_Q)+0.2)
    lowQ_ = 10**(np.log10(med1_Q)-0.2)
    
    ax2[1].plot(x1_SF,upSF_,color='dodgerblue',ls=':',lw=4)
    ax2[1].plot(x1_SF,lowSF_,color='dodgerblue',ls=':',lw=4)  
    ax2[1].plot(x1_Q,upQ_,color='orangered',ls=':',lw=4)
    ax2[1].plot(x1_Q,lowQ_,color='orangered',ls=':',lw=4)
    
    
    
    ax2[2].plot(x2_SF,med2_SF, ls='-',  color='dodgerblue',lw=5)
    ax2[2].fill_between(x2_SF,up2_SF,down2_SF, color='dodgerblue',alpha=0.4)
    ax2[2].plot(x2_Q,med2_Q, ls='-',  color='orangered',lw=5)
    ax2[2].fill_between(x2_Q,up2_Q,down2_Q, color='orangered',alpha=0.4)    
    
    upSF_ = 10**(np.log10(med2_SF)+0.2)
    lowSF_ = 10**(np.log10(med2_SF)-0.2)
    upQ_ = 10**(np.log10(med2_Q)+0.2)
    lowQ_ = 10**(np.log10(med2_Q)-0.2)
    
    ax2[2].plot(x2_SF,upSF_,color='dodgerblue',ls=':',lw=4)
    ax2[2].plot(x2_SF,lowSF_,color='dodgerblue',ls=':',lw=4)  
    ax2[2].plot(x2_Q,upQ_,color='orangered',ls=':',lw=4)
    ax2[2].plot(x2_Q,lowQ_,color='orangered',ls=':',lw=4)
'''

    ############ BT selected ###################

    fig,(ax1,ax2,ax3) = plt.subplots(3,3,figsize = (32,32),sharey=True,sharex=True)

    z=1
    x0_SF,med0_SF,up0_SF,down0_SF = np.loadtxt('./SAM/bt_selected/ltgs/z'+str(z)+'/median.dat', unpack=True)
    x0_Q,med0_Q,up0_Q,down0_Q = np.loadtxt('./SAM/bt_selected/etgs/z'+str(z)+'/median.dat', unpack=True)

   # m0 = np.ma.masked_greater(x0_Q,50).mask

    #x0_Q = x0_Q[m0]
    
    #med0_Q,up0_Q,down0_Q = make_savgol(med0_Q[m0],up0_Q[m0],down0_Q[m0], window=3)

    z=3
    x1_SF,med1_SF,up1_SF,down1_SF = np.loadtxt('./SAM/bt_selected/ltgs/z'+str(z)+'/median.dat', unpack=True)
    x1_Q,med1_Q,up1_Q,down1_Q = np.loadtxt('./SAM/bt_selected/etgs/z'+str(z)+'/median.dat', unpack=True)

    #m1 = np.ma.masked_greater(x1_Q,50).mask

    #x1_Q = x1_Q[m1]
   # med1_Q,up1_Q,down1_Q = make_savgol(med1_Q[m1],up1_Q[m1],down1_Q[m1],window=5)
    
    
    z=5
    x2_SF,med2_SF,up2_SF,down2_SF = np.loadtxt('./SAM/bt_selected/ltgs/z'+str(z)+'/median.dat', unpack=True)
    x2_Q,med2_Q,up2_Q,down2_Q = np.loadtxt('./SAM/bt_selected/etgs/z'+str(z)+'/median.dat', unpack=True)

    #m2 = np.ma.masked_greater(x2_Q,50).mask
    #x2_Q = x2_Q[m2]
    #med2_Q = med2_Q[m2]
    #up2_Q = up2_Q[m2]
    #down2_Q = down2_Q[m2]



    ax3[0].plot(x0_SF,med0_SF, ls='-',  color='dodgerblue',lw=5)
    ax3[0].fill_between(x0_SF,up0_SF,down0_SF, color='dodgerblue',alpha=0.4)
    ax3[0].plot(x0_Q,med0_Q, ls='-',  color='orangered',lw=5)
    ax3[0].fill_between(x0_Q,up0_Q,down0_Q, color='orangered',alpha=0.4)
    
    upSF_ = med0_SF +0.2# 10**(np.log10(med0_SF)+0.2)
    lowSF_ = med0_SF-0.2#10**(np.log10(med0_SF)-0.2)
    upQ_ = med0_Q+0.2
    lowQ_ = med0_Q-0.2
    
    ax3[0].plot(x0_SF,upSF_,color='dodgerblue',ls=':',lw=4)
    ax3[0].plot(x0_SF,lowSF_,color='dodgerblue',ls=':',lw=4)    
    ax3[0].plot(x0_Q,upQ_,color='orangered',ls=':',lw=4)
    ax3[0].plot(x0_Q,lowQ_,color='orangered',ls=':',lw=4)
    
    
    ax3[1].plot(x1_SF,med1_SF, ls='-',  color='dodgerblue',lw=5)
    ax3[1].fill_between(x1_SF,up1_SF,down1_SF, color='dodgerblue',alpha=0.4)
    ax3[1].plot(x1_Q,med1_Q, ls='-',  color='orangered',lw=5)
    ax3[1].fill_between(x1_Q,up1_Q,down1_Q, color='orangered',alpha=0.4)    

    
    upSF_ = med1_SF+0.2
    lowSF_ = med1_SF-0.2
    upQ_ = med1_Q+0.2
    lowQ_ = med1_Q-0.2
    
    ax3[1].plot(x1_SF,upSF_,color='dodgerblue',ls=':',lw=4)
    ax3[1].plot(x1_SF,lowSF_,color='dodgerblue',ls=':',lw=4)  
    ax3[1].plot(x1_Q,upQ_,color='orangered',ls=':',lw=4)
    ax3[1].plot(x1_Q,lowQ_,color='orangered',ls=':',lw=4)
    
    
    
    ax3[2].plot(x2_SF,med2_SF, ls='-',  color='dodgerblue',lw=5)
    ax3[2].fill_between(x2_SF,up2_SF,down2_SF, color='dodgerblue',alpha=0.4)
    ax3[2].plot(x2_Q,med2_Q, ls='-',  color='orangered',lw=5)
    ax3[2].fill_between(x2_Q,up2_Q,down2_Q, color='orangered',alpha=0.4)    
    
    upSF_ = med2_SF+0.2
    lowSF_ = med2_SF-0.2
    upQ_ = med2_Q+0.2
    lowQ_ = med2_Q-0.2
    
    ax3[2].plot(x2_SF,upSF_,color='dodgerblue',ls=':',lw=4)
    ax3[2].plot(x2_SF,lowSF_,color='dodgerblue',ls=':',lw=4)  
    ax3[2].plot(x2_Q,upQ_,color='orangered',ls=':',lw=4)
    ax3[2].plot(x2_Q,lowQ_,color='orangered',ls=':',lw=4)





    xx = np.log10(np.arange(50,500))
    for i in range(3):
        ax1[i].plot(xx,np.log10(0.015) + xx, color='black', lw=2,ls='--', label='Kravtsov 2013')
        ax2[i].plot(xx,np.log10(0.015)+xx,color='black', ls='--',lw=2)
        ax3[i].plot(xx,np.log10(0.015)+xx,color='black', ls='--',lw=2)
    
    
  #  for i in range(0,3):
  #      ax1[i].set_xscale('log')
  #      ax2[i].set_xscale('log')
  #      ax3[i].set_xscale('log')

   #     ax1[i].set_yscale('log')
   #     ax2[i].set_yscale('log')
   #     ax3[i].set_yscale('log')
 
    
    
    
   # ax1[0].get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
   # ax1[0].minorticks_off()
   # ax1[0].set_yticks([2,4,10,20,40])

    #ax2[0].get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    #ax2[0].minorticks_off()
    #ax2[0].set_yticks([2,4,10,20,40])

    #ax3[0].get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    #ax3[0].minorticks_off()
    #ax3[0].set_yticks([2,4,10,20,40])
    
#    ax1[1].get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
#    ax1[1].minorticks_off()
#    ax1[1].set_yticks([2,4,10,20,40])
    
#    ax1[2].get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
#    ax1[2].minorticks_off()
#    ax1[2].set_yticks([2,4,10,20,40])

    #for i in range(0,3):
    #    
    #    ax3[i].get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    #    ax3[i].minorticks_off()
    #    ax3[i].set_xticks([50,100,300])

    #    ax3[0].set_xlim(50)
    #    ax3[1].set_xlim(50)
    #    ax3[2].set_xlim(50)

    #ax1[0].set_ylim(0.5)
    #ax2[0].set_ylim(0.5)
    #ax3[0].set_ylim(0.5)
    
    ax1[0].set_title('z=0')
    ax1[1].set_title('z=1')
    ax1[2].set_title('z=2')
    
    ax1[0].set_ylabel('$R_e \ [kpc] $')
    ax2[0].set_ylabel('$R_e \ [kpc] $')
    ax3[0].set_ylabel('$R_e \ [kpc] $')
    
    ax3[0].set_xlabel('$R_h \ [kpc]$')
    ax3[1].set_xlabel('$R_h \ [kpc]$')
    ax3[2].set_xlabel('$R_h \ [kpc]$')
    
    ax1[0].legend( frameon=False, fontsize=35)
    plt.subplots_adjust(wspace=0, left=0.1, hspace=0,right=0.8, )

    
    fig.text(0.82,0.75,'Illustris TNG', fontsize=50)
    fig.text(0.82,0.5,'Rome SAM \n sSFR selected', fontsize=50)
    fig.text(0.82,0.25,'Rome SAM \n B/T selected', fontsize=50)


    #plt.tight_layout()
    plt.savefig('/home/lz1f17/Pictures/Paper/TNGsizes_new.pdf')
       
       
