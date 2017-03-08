import emcee
import matplotlib as mpl
mpl.use('TkAgg')
import seaborn as sea
import matplotlib.pyplot as plt
import corner
import scipy.stats as stats
import numpy as np
import astropy.units as u
import astropy.wcs as wcs
import astropy.io.fits as fits
import datetime
import os
import math
from IPython import display


default_crpix1=[5236.283490922594,3162.773007622719,1051.069038778534,-1062.272650429601,-3157.2377640915056,5236.1366766101019,\
    3163.9128920108824,1053.6999541273758,-1069.7269810065368,-3158.5262939764566]

default_crpix2=[-118.8920592518456,-118.7585232944771,-119.0352833407502,-116.8035616903668,-117.16262318037234,4085.7234028874109,\
    4101.2911391270209,4107.5698514917149,4099.2096669553466,4087.898125411612]

default_crota2=[0.1903798950331974,0.2082959367048585,0.1971478000328974,0.2098078872741506,\
    0.22713412160210628,0.22198114823159171,0.22452571939824689,0.23341475817443058,0.2333373458198231,0.20265102927000198]

default_cdelt1 = -5.609310e-05
default_cdelt2 = 5.609582e-05

class SuprimeCamMCMC(object):
    """
    Class to perform Bayesian inference for SuprimeCam astrometric parameters.

    Parameters:
        data: i x 4 x k numpy array, with i sources for each of k chips.
            Four rows are RA, DEC, X, and Y.

        image_filename: .fits file to take WCS information for crval1, crval2, and angle.
        Should be tmos_ file for specific object for standard scampipe naming conventions.

        crval1, crval2, angle: initial parameters if .fits file not provided

        cdelt1, cdelt2: initial values for scale parameters

        crpix1,crpix2,crota2: initial values for the chip-specific quantities

        single_chip: if n_chips = 1, specify which chip we are fitting

        mix: if mix specified, treat model as mixture of outliers and real sources

        init_good_frac: initial value for good source fraction.

        diag_plots_dir: directory to output diagnostic plots to

        image_dir: directory where images (and catalogs for the moment) live
    """
    def __init__(self, data, n_walkers = 100, n_chips=1,\
        init_crpix1=default_crpix1,\
        init_crpix2=default_crpix2,\
        init_crota2=default_crota2,\
        init_cdelt1=default_cdelt1,\
        init_cdelt2=default_cdelt2,\
        init_crval1=0.0,init_crval2=0.0,init_angle=0.0,\
        image_filename=None,\
        fix_crvals=True,fix_cdelts=False,fix_angle=True,
        x_pix_sig=1.,y_pix_sig=1.,single_chip=1,
        x_out_pix_sig=5,y_out_pix_sig=5,
        mix=True,init_good_frac=0.95,diag_plots=True,
        diag_plots_dir='mcmc_diag_plots/',
        image_dir=''):

        self.data = data
        self.n_walkers = n_walkers
        self.n_chips = n_chips
        self.x_pix_sig = x_pix_sig
        self.y_pix_sig = y_pix_sig
        self.image_dir = image_dir

        #we might want to use a mixture model to account for outliers.
        #if so, we record the sigmas to use for the outlier distribution
        self.x_out_pix_sig = x_out_pix_sig
        self.y_out_pix_sig = y_out_pix_sig
        self.mix = mix

        if self.n_chips < 2:
            self.chip = single_chip

        #infer the values for globals from a .fits file, if provided
        if image_filename is not None:
            init_crval1,init_crval2,init_cdelt1,init_cdelt2,init_crota2 = \
                self.get_wcs(image_filename)

        #check if we're fixing the global parameters and save values if true
        self.fix_crvals = fix_crvals
        if self.fix_crvals:
            self.crval1 = init_crval1
            self.crval2 = init_crval2

        self.fix_cdelts = fix_cdelts
        if self.fix_cdelts:
            self.cdelt1 = init_cdelt1
            self.cdelt2 = init_cdelt2

        self.fix_angle = fix_angle
        if self.fix_angle:
            self.angle = init_angle

        #pack the starting values into the intial theta
        self.theta_init = self.theta_pack(init_crpix1,init_crpix2,init_crota2,init_cdelt1,\
            init_cdelt2,init_crval1,init_crval2,init_angle,init_good_frac)

        print('Initial WCS Image Values: ',self.theta_init)

        #initialize the WCS objects

        self.init_wcs(image_filename)
        #print(image_filename)
        #find the dimensionality of the problem
        self.n_dim = 3*self.n_chips + 2*int(not self.fix_crvals) + 2*int(not self.fix_cdelts) \
            + int(not self.fix_angle) + int(self.mix)

        #initialize the sampler
        self.sampler = emcee.EnsembleSampler(self.n_walkers,self.n_dim,self.ln_like)

        #shake up the thetas to initialize the walkers
        self.pos = np.array([(self.theta_init \
            + np.array(self.theta_init)*1.e-4*np.random.randn(self.n_dim))\
             for j in range(self.n_walkers)])

        #check if making diagnostic plots, and make directory if it doesn't exist already
        self.diag_plots = diag_plots
        if self.diag_plots:
            self.diag_plots_dir = diag_plots_dir
            if not os.path.isdir(self.diag_plots_dir):
                os.makedirs(self.diag_plots_dir)


    def get_wcs(self,filename):
        '''
        Get the global values from a specified fits file.
        '''
        hdu_list = fits.open(self.image_dir+filename)
        crval1 = hdu_list[0].header['crval1']
        crval2 = hdu_list[0].header['crval2']
        cdelt1 = hdu_list[0].header['cdelt1']
        cdelt2 = hdu_list[0].header['cdelt2']
        try:
            angle = np.array([math.degrees(math.acos(hdu_list[0].header['cd1_1'] / hdu_list[0].header['cdelt1']))])

        #Sometimes this produces a quantity greater than one. This should only happen when angle is ~0.
        except ValueError:
            print('Error reading initial angle, setting initial angle to zero.')
            angle = np.array([0.00001])
        except KeyError:
        #If WCS is not in matrix form, crota2 header should exist, so take this instead.
            angle = np.array([hdu_list[0].header['crota2']])
        #angle=0.0

        return (crval1,crval2,cdelt1,cdelt2,angle)

    def theta_pack(self,crpix1,crpix2,crota2,cdelt1,cdelt2,crval1,crval2,angle,good_frac):
        """
        Wrap parameters into 1D array for emcee.
        """
        theta = []

        #for fitting all chips, loop over all chips
        if self.n_chips < 1:
            for i in range(self.n_chips):
                theta.append(crpix1[i])

            for i in range(self.n_chips):
                theta.append(crpix2[i])

            for i in range(self.n_chips):
                theta.append(crota2[i])

        else:
            theta.append(crpix1[self.chip])
            theta.append(crpix2[self.chip])
            theta.append(crota2[0])


        #check values we might want to fix, and add them if they're left free
        if not self.fix_cdelts:
            theta.append(cdelt1)
            theta.append(cdelt2)

        if not self.fix_crvals:
            theta.append(crval1)
            theta.append(crval2)

        if not self.fix_angle:
            theta.append(angle)

        #if we're mixing, need to include the outlier probability
        if self.mix:
            theta.append(good_frac)

        return theta

    def theta_unpack(self,theta):
        """
        Unpack the 1D theta array for easier handling of parameters.

        Sneaky: if global parameters are being fixed, we slip the
        fixed values in here.
        """
        crpix1 = theta[:self.n_chips]
        crpix2 = theta[self.n_chips:2*self.n_chips]
        crota2 = theta[2*self.n_chips:3*self.n_chips]

        ind = 3*self.n_chips
        if not self.fix_cdelts:
            cdelt1 = theta[ind]
            cdelt2 = theta[ind+1]
            ind += 2
        else:
            cdelt1 = self.cdelt1
            cdelt2 = self.cdelt2

        if not self.fix_crvals:
            crval1 = theta[ind]
            crval2 = theta[ind+1]
            ind += 2
        else:
            crval1 = self.crval1
            crval2 = self.crval2

        if not self.fix_angle:
            angle = theta[ind]
            ind += 1
        else:
            angle = self.angle

        if self.mix:
            good_frac = theta[ind]
            ind += 1
        else:
            good_frac = 1.

        return (crpix1,crpix2,crota2,cdelt1,cdelt2,crval1,crval2,angle,good_frac)

    def init_wcs(self,filename):
        """initialize the WCS objects for each chip"""
        self.wcs = []
        for i in range(self.n_chips):
            #self.wcs.append(wcs.WCS(naxis=2))
            #self.wcs[i].WCS.ctype= ["RA---TAN", "DEC--TAN"]
            hdu_list = fits.open(self.image_dir+filename)
            hdu_list[0].header.pop('cd1_1',None)
            hdu_list[0].header.pop('cd1_2',None)
            hdu_list[0].header.pop('cd2_1',None)
            hdu_list[0].header.pop('cd2_2',None)
            self.wcs.append(wcs.WCS(hdu_list[0].header))
            #del self.wcs[i].wcs.cd

    def ln_like(self,theta):
        """
        Returns the log of the likelihood for parameters theta.
        """
        crpix1,crpix2,crota2,cdelt1,cdelt2,crval1,crval2,angle,good_frac\
            = self.theta_unpack(theta)

        self.x = []
        self.y = []

        #loop over each chip and perform the transforms from ra/dec to x,y
        for i in range(self.n_chips):
            self.wcs[i].wcs.crpix = [crpix1[i],crpix2[i]]
            self.wcs[i].wcs.crval = [crval1,crval2]
            self.wcs[i].wcs.cdelt = np.array([cdelt1,cdelt2])

            #this rotation angle isn't correct (I think it's this + angle?)
            self.wcs[i].wcs.crota = [0.,crota2[i] + angle]

            #transform RA, DEC to appropriate X, Y.
            self.new_xy = self.wcs[i].wcs_world2pix(self.data[:,2:4,i],1)
            self.x.append(self.new_xy[:,0])
            self.y.append(self.new_xy[:,1])

        self.x_diffs = self.data[:,0,:].ravel() - self.x
        self.y_diffs = self.data[:,1,:].ravel() - self.y

        #now that differences are calculated, calculate logPDFs
        self.log_like_x = np.log(good_frac)+\
            stats.norm.logpdf(self.x_diffs,loc=0.,scale=self.x_pix_sig)
        self.log_like_y = np.log(good_frac)+\
            stats.norm.logpdf(self.y_diffs,loc=0.,scale=self.y_pix_sig)

        #Calculate outlier logPDFs (which will be zero if )
        self.log_like_out_x = (np.log(1.-good_frac))+\
            stats.norm.logpdf(self.x_diffs,loc=0.,scale=self.x_out_pix_sig)
        self.log_like_out_y = (np.log(1.-good_frac))+\
            stats.norm.logpdf(self.y_diffs,loc=0.,scale=self.y_out_pix_sig)

        #Calculate prior for our specific theta
        ln_prior = self.ln_prior(theta)
        if not np.isfinite(ln_prior):
            return -np.inf

        #add the outlier PDFs and the
        self.full_like_x = np.sum(np.logaddexp(self.log_like_x,self.log_like_out_x))
        self.full_like_y = np.sum(np.logaddexp(self.log_like_y,self.log_like_out_y))


        log_like = self.full_like_x + self.full_like_y+ ln_prior

        return log_like

    def ln_prior(self,theta):
        """calculate the log of the prior"""
        crpix1,crpix2,crota2,cdelt1,cdelt2,crval1,crval2,angle,good_frac\
         = self.theta_unpack(theta)

        crpix1_i,crpix2_i,crota2_i,cdelt1_i,cdelt2_i,crval1_i,crval2_i,angle_i,good_frac_i =\
             self.theta_unpack(self.theta_init)

        if np.any(np.abs(np.array(crpix1) - np.array(crpix1_i)) > 1000):
            return -np.inf
        if np.any(np.abs(np.array(crpix2) - np.array(crpix2_i)) > 1000):
            return -np.inf
        if np.any(np.abs(np.array(crota2) - np.array(crota2_i)) > 5.):
            return -np.inf

        if (-6.0e-5 > cdelt1) or (cdelt1 > -5.0e-5):
            return -np.inf

        if (6.0e-5 < cdelt2) or (cdelt2 < 5.0e-5):
            return -np.inf

        if self.mix:
            mix_prior = stats.beta.logpdf(good_frac,0.8,0.2)
        else:
            mix_prior = 0.

        return mix_prior


    def run_sampler(self,nsteps,progress_bar=False):
        """take n mcmc steps"""
        if progress_bar:
            for i, (self.pos, lnp, state) in enumerate(self.sampler.sample(self.pos, iterations=nsteps)):
                if (i+1) % 100 == 0:
                    print("{0:.1f}%".format(100 * float(i) / nsteps))

        else:
            self.sampler.run_mcmc(self.pos,nsteps)

        self.chain = self.sampler.chain
        self.pos = self.sampler.chain[:,-1,:]


    def get_stats(self,nburn=0,hi_region=16.,lo_region=84.,print_stats=True):
        """
        Print and return descriptive statistics for the walker chains

        hi_region: precentile (0-100) for high credible region estiamte
        lo_region: precentile (0-100) for low credible region estimate
        """
        meds = []
        lo = []
        hi = []
        for i in range(self.n_dim):
            chain = self.sampler.chain[:,nburn:,i].flatten()
            meds.append(np.median(chain))
            lo.append(np.percentile(chain,lo_region))
            hi.append(np.percentile(chain,hi_region))

            if print_stats:
                print 'Statistics for parameter '+str(i)+'.'
                print 'Median: '+str(meds[i])+', Low: '+str(lo[i])+\
                    ', High: '+str(hi[i])

        self.meds = meds

        return (meds,lo,hi)

    def donate_astrometry(self,filename,meds=None,pre='aS',clobber=False):
        """
        Donate updated astrometric solution to the provided filename.

        Meds: Array-like, consisting of [crpix1,crpix2,crota2,crdelt1,crdelt2]
        If not provided, will attempt to use self.meds (created if get_stats has been run)

        Pre: prefix for the updated astrometry file. Default is consistent for scampipe

        Clobber: if true, clobber donated astrometry files
        """
        hdu = fits.open(self.image_dir+filename)

        header = hdu[0].header
        if meds is None:
            meds = self.meds

        header['crpix1'] = meds[0]
        header['crpix2'] = meds[1]
        header['crota2'] = meds[2]
        header['cdelt1'] = meds[3]
        header['cdelt2'] = meds[4]
        header['crval1'] = self.wcs[0].wcs.crval[0]
        header['crval2'] = self.wcs[0].wcs.crval[1]

		#take out the cd keywords which will confuse the WCS solution otherwise
        header.pop('cd1_1',None)
        header.pop('cd1_2',None)
        header.pop('cd2_1',None)
        header.pop('cd2_2',None)
		
		#in old data, the keywords are stored in a different format, take those out too.
        header.pop('PC001001',None)
        header.pop('PC001002',None)
        header.pop('PC002001',None)
        header.pop('PC002002',None)

        header['History'] = 'Updated Astrometry Keywords :' + str(datetime.datetime.now())

        hdu[0].header = header

        hdu.writeto(self.image_dir+pre+filename,overwrite=clobber)

    def make_diag_plots(self,nburn=0,trace_burn=False,trace_walkers=30,file_root=''):
        """
        Make trace, corner, and residual plots for the chain.
        """
        display.clear_output(wait=True)
        display.display(plt.gcf())

        fig,ax = plt.subplots(figsize=(17,8),ncols=3,nrows=2,dpi=1000)

        if trace_burn:
            nburn_trace = nburn
        else:
            nburn_trace = 0

        for i in range(self.n_dim):
            for j in range(trace_walkers):
                ax[i/3,i%3].plot(self.sampler.chain[j,nburn_trace:,i])


        fig.savefig(self.diag_plots_dir+file_root+'_trace_plots'+'.pdf')
        plt.close(fig)

        fig = corner.corner(self.sampler.chain[:,nburn:,:].reshape((-1,self.n_dim)))
        fig.savefig(self.diag_plots_dir+file_root+'_corner_plots'+'.pdf')
        plt.close(fig)

        if hasattr(self, 'meds'):
            meds = self.meds
        else:
            meds,lows,highs = self.get_stats(nburn=nburn,print_stats=False)

        new_wcs = self.wcs[0]

        new_wcs.crota = [0.,meds[2]]
        new_wcs.crpix1 = meds[0]
        new_wcs.crpix2 = meds[1]

        new_xy = self.wcs[0].wcs_world2pix(self.data[:,2:4,0],1)
        new_x = new_xy[:,0]
        new_y = new_xy[:,1]
        x_diffs = new_x - self.data[:,0,0]
        y_diffs = new_y - self.data[:,1,0]

        fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(12,5))
        ax[0,0].scatter(self.data[:,0],x_diffs)
        ax[0,1].scatter(self.data[:,1],y_diffs)

        ax[1,0].scatter(self.data[:,1],x_diffs)
        ax[1,1].scatter(self.data[:,0],y_diffs)

        fig.savefig(self.diag_plots_dir+file_root+'_residual_plots'+'.pdf')
        plt.close(fig)


















