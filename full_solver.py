scam_chip_names = ['chihiro','clarisse','fio','kiki','nausicaa','ponyo','san','satsuki','sheeta','sophie']
scam_chip_names_old = ['w67c1', 'w6c1', 'si006s', 'si002s', 'w7c3', 'w93c2', 'w9c2', 'si005s', 'si001s', 'w4c5']

from suprimecam_mcmc import SuprimeCamMCMC
import numpy as np
import sys

class SCamMultiChipSolver(object):
    """
    Class to perform several MCMC iterations for SCam chips.

    filelist: list of image files to loop over. Must follow standard
    scampipe formatting.
    """
    def __init__(self, filelist, image_dir=''):
        self.filelist = np.genfromtxt(image_dir+filelist,dtype=str)
        self.mcmcs = []
        self.meds= {}
        self.image_dir = image_dir


    def fit_images(self,filelist=None,nsteps=800,nburn=500,make_diag_plots=True,\
        start_chip=0,end_chip=-1,**kwargs):
        """
        Fit astrometry for a file or list of files.
        nsteps: number of steps for each sampler to take
        nburn: number of steps to burn for inference of parameters

        make_diag_plots: if true, output plots to mcmc_diag_plots/
        currently pretty glitchy, will also output plots to standard viewer
        and introduces some overhead

        start_chip: index of first file to loop over
        end_chip: index of last file to loop over (-1 will loop through last file)

        """
        if filelist is None:
            filelist = self.filelist


        #loop over the list of files
        for file in filelist[start_chip:end_chip]:
            #split the filename
            split = file.split("_",2)
            image_designation = split[1]

            chipname = split[2].split(".",2)[0]
            try:
                chipnumber = scam_chip_names.index(chipname)
            except ValueError:
                chipnumber = scam_chip_names_old.index(chipname)

            cat_name = image_designation+'_'+chipname+'.csv'
            #cat_name = 'sdss_matches_chip_'+str(chipnumber)+'.cat'

            #read in the data
            data = self.read_data(cat_name)

            tmos = "tmos_" + image_designation + ".fits"
            print tmos
            print 'Solving astrometry for image '+image_designation+', chip '+chipname

            this_mcmc = SuprimeCamMCMC(data,image_filename=tmos,single_chip=chipnumber,
                image_dir=self.image_dir,**kwargs)
            this_mcmc.run_sampler(nsteps)
            meds,lows,highs = this_mcmc.get_stats(nburn=nburn)

            this_mcmc.donate_astrometry(file,clobber=True)

            if make_diag_plots:
                this_mcmc.make_diag_plots(nburn=nburn,file_root=image_designation+'_'+chipname)

            self.meds[image_designation+chipname] = meds

            self.mcmcs.append(this_mcmc)

            print 'Solve finished, diagnostic plots output'

    def read_data(self,filename):
        """
        Function to read in data from a catalog file.

        Assumes formatting to be x,y,ra,dec
        """
        data = np.genfromtxt(self.image_dir+filename,delimiter=',')
        data = data.reshape(data.shape[0],data.shape[1],1)

        return data


def main():
    full = SCamMultiChipSolver(sys.argv[1])
    full.fit_images(nsteps=600,nburn=400,make_diag_plots=True,start_chip=0,end_chip=None)

if __name__ == "__main__":
    main()















