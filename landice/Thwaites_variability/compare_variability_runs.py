#!/usr/bin/env python
'''
Script to compare some scalar values from different runs of Thwaites melt variability experiment.
'''

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm

outfname = 'globalStats.nc'
runs=[ adir for adir in sorted(os.listdir('.')) if (os.path.isdir(adir) and os.path.isfile(os.path.join(adir, outfname)))]
print "Original run list:", runs


# get just the 0 phase for each amp/per
#runs[:] = [r for r in runs if 'pha0.00' in r]
# get all the phase for one amp/per
#runs[:] = [r for r in runs if 'amp300_per70' in r]
#runs.append('steady')

# reorder to put the 'control' runs at the beginning
special_runs = ('steady', 'no-melt')
for r in special_runs:
   if r in runs:
      runs.remove(r)
      runs.insert(0, r)

# optionally exclude some subset
#runs[:] = [r for r in runs if not 'amp300_per02' in r]

print "Will process the following directories: ", runs

# --- Define some needed parameters
rhoi = 910.0  # ice density
windowLength = 7      # window length for smoothing some time series
# ----------------

# ------ Needed functions ------

def xtime2numtime(xtime):
  """Define a function to convert xtime character array to numeric time values using datetime objects"""
  # First parse the xtime character array into a string
  xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function

  dt = []
  for stritem in xtimestr:
      itemarray = stritem.strip().replace('_', '-').replace(':', '-').split('-')  # Get an array of strings that are Y,M,D,h,m,s
      results = [int(i) for i in itemarray]
      if (results[0] < 1900):  # datetime has a bug where years less than 1900 are invalid on some systems
         results[0] += 1900
      dt.append( datetime.datetime(*results) ) # * notation passes in the array as arguments

  numtime = netCDF4.date2num(dt, units='seconds since '+str(dt[0]))   # use the netCDF4 module's function for converting a datetime to a time number
  return numtime / (3600.0*24.0*365.0) # in years

def xtime2numtimeMy(xtime):
  """Define a function to convert xtime character array to numeric time values using local arithmetic"""
  # First parse the xtime character array into a string
  xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function

  numtime = np.zeros( (len(xtimestr),) )
  ii = 0
  dayOfMonthStart = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
  for stritem in xtimestr:
      itemarray = stritem.strip().replace('_', '-').replace(':', '-').split('-')  # Get an array of strings that are Y,M,D,h,m,s
      results = [int(i) for i in itemarray]
      numtime[ii] = results[0] + (dayOfMonthStart[results[1]-1]-1 + results[2]) / 365.0  # decimal year
      ii += 1
  return numtime

def GTtoSL(GT):
   #return (GT-VAF[0])*1.0e12/1028.0/1000.0**3/362.0e6*1000.0*1000.0
   return GT *1.0e12/1028.0/1000.0**3/362.0e6*1000.0*1000.0 * -1.0

# --- Define data structures ---

class modelRun:
   def __init__(self, run):
      '''
      This reads results from a model run and saves and analyzes the needed results.

      run = name of subdir in which the run was performed
      '''

      # some metadata about run
      if 'amp' in run:
         self.amp = float(run[3:6])
         self.per = float(run[10:12])
         self.pha = float(run[16:])
      else:
         self.amp = 0.0
         self.per = 0.0
         self.pha = 0.0

      #print self.amp, self.per, self.pha

      f = netCDF4.Dataset(run + '/' + outfname, 'r')
      self.nt = len(f.dimensions['Time'])
      self.yrs = np.zeros((self.nt,))
      #yrs = f.variables['daysSinceStart'][:] / 365.0
      self.xtimes = f.variables['xtime'][:]
   
   
      self.yrs = xtime2numtimeMy(self.xtimes)
      self.dyrs = self.yrs[1:] - self.yrs[0:-1]
   

      # resampled version of time array - needed for filtering.
      # (Filtering needed b/c the occasional tiny time step that the model is exhibiting leads to inaccurate (noisy) derivatives)
      resampEndtime = self.yrs.max()
      resampEndtime = 500.0
      self.resampYrs = np.linspace(0.0, resampEndtime, num=resampEndtime*12*2)
      self.dresampYrs = self.resampYrs[1:] - self.resampYrs[0:-1]
   
      self.melt = f.variables['avgSubshelfMelt'][:]
      self.resampMelt = np.interp(self.resampYrs, self.yrs, self.melt) # generate y values for each x
      self.totalmelt = f.variables['totalFloatingBasalMassBal'][:] / 1.0e12
      self.cumumelt = self.totalmelt*0.0
      for t in np.arange(self.nt-1)+1:
         self.cumumelt[t] = self.cumumelt[t-1] + self.totalmelt[t-1] * self.dyrs[t-1]
   
      self.VAF = f.variables['volumeAboveFloatation'][:] / 1.0e12 * rhoi
      self.resampVAF = np.interp(self.resampYrs, self.yrs, self.VAF) # generate y values for each x
      self.VAFsmooth = scipy.signal.savgol_filter(self.resampVAF, window_length=windowLength, polyorder=1)
      #VAFrate = (VAF[1:]-VAF[0:-1] ) / dyrs
      self.VAFsmoothrate = (self.VAFsmooth[1:] - self.VAFsmooth[0:-1] ) / self.dresampYrs
   
      self.GA = f.variables['groundedIceArea'][:] / 1.0e3**2
      self.resampGA = np.interp(self.resampYrs, self.yrs, self.GA) # generate y values for each x
      self.GAsmooth = scipy.signal.savgol_filter(self.resampGA, window_length=windowLength, polyorder=1)
      self.GArate = (self.GAsmooth[1:] - self.GAsmooth[0:-1]) / self.dresampYrs

      # calculate lag
      #self.VAFeven = np.arange(self.VAF[0], 100000.0, -250.0)
      #self.VAFeven = np.linspace(70000.0, self.VAF[0], 350)
      self.VAFeven = np.linspace(45000.0, 209053.0, 2000)  #229728878784122
      self.timeOnVAFeven = np.interp(self.VAFeven, self.VAF[::-1], self.yrs[::-1])
      self.VAFrateOnVAFeven = np.interp(self.VAFeven, self.VAFsmooth[:0:-1], self.VAFsmoothrate[::-1])
      self.meltOnVAFeven = np.interp(self.VAFeven, self.VAF[::-1], self.melt[::-1])
      # redo as SLR
      #self.VAFeven = GTtoSL(np.linspace(70000.0, 209053.0, 2000) - 209053.0)  #229728878784122
      #self.timeOnVAFeven = np.interp(self.VAFeven, GTtoSL(self.VAF[::-1] - self.VAF[0]), self.yrs[::-1])
      #self.VAFrateOnVAFeven = np.interp(self.VAFeven, self.VAFsmooth[:0:-1], self.VAFsmoothrate[::-1])
      #self.meltOnVAFeven = np.interp(self.VAFeven, self.VAF[::-1], self.melt[::-1])

      self.GAeven = np.linspace(79737808225.442/1.0e3**2, 173169904261.254/1.0e3**2, 500)
      self.timeOnGAeven = np.interp(self.GAeven, self.GA[::-1], self.yrs[::-1])
      #self.VAFrateOnGAeven = np.interp(self.GAeven, self.VAFsmooth[:0:-1], self.VAFsmoothrate[::-1])
      self.meltOnGAeven = np.interp(self.GAeven, self.GA[::-1], self.melt[::-1])


      self.GLflux = f.variables['groundingLineFlux'][:] / 1.0e12

      self.floatingVol = f.variables['floatingIceVolume'][:]# / 1.0e12 * rhoi
      self.floatingArea = f.variables['floatingIceArea'][:]

      f.close()


# ================
# Loop over runs and collect needed data
# ================
runData = {}  # init empty dictionary
groups = {} # groups are a dictionary of dictionaries.  May want to changes this to a dictionary of custom objects that can include metadata about the group alongside the group's data children.  KISS for now.
for run in runs:
   print "Processing run: " + run

   # Build a list that contains all run data objects
   runData[run] = modelRun(run)

   # Populate each run into its group for ensemble runs (not for controls)
   if 'amp' in run:
      groupName = run[0:12]  # this is the amp and per
      #print groupName
      phase = run[13:]
      #print phase
      if groupName not in groups:
         groups[groupName] = {} # init an empty dict to add group members
      groups[groupName][phase] = runData[run]  # stick the run object into this group

print "Processing complete.\n"

#print groups


# --- set up figure axes ---

print "Setting up figure axes."

fig = plt.figure(1, facecolor='w')

nrow=4
ncol=2

xtickSpacing = 50.0

# melt forcing
axMeanMelt = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('Year')
plt.ylabel('mean melt (m/yr)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

axTotalMelt = fig.add_subplot(nrow, ncol, 3, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('total melt (Gt/yr)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

axCumuMelt = fig.add_subplot(nrow, ncol, 5, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('cumulative melt (Gt)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

# VAF
axVAF = fig.add_subplot(nrow, ncol, 2, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF (Gt)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

axVAFrate = fig.add_subplot(nrow, ncol, 4, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF rate (Gt/yr)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

# grounded area
axGA = fig.add_subplot(nrow, ncol, 6, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('Grounded area (km^2)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

axGArate = fig.add_subplot(nrow, ncol, 8, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('Grounded area rate (km^2/yr)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()



# =================== repeat cleaner version for presentation ===========

fig2 = plt.figure(2, facecolor='w')

nrow=5
ncol=1

# melt forcing
ax2MeanMelt = fig2.add_subplot(nrow, ncol, 1)
plt.xlabel('Year')
plt.ylabel('mean melt (m/yr)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

# VAF
ax2VAF = fig2.add_subplot(nrow, ncol, 2, sharex=ax2MeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF (Gt)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

ax2VAFrate = fig2.add_subplot(nrow, ncol, 3, sharex=ax2MeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF rate (Gt/yr)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

ax2GLflux = fig2.add_subplot(nrow, ncol, 4, sharex=ax2MeanMelt)
plt.xlabel('Year')
plt.ylabel('GL flux (Gt/yr)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

ax2floatVol = fig2.add_subplot(nrow, ncol, 5, sharex=ax2MeanMelt)
plt.xlabel('Year')
plt.ylabel('mean shelf thickness (m)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()


# ======
# this figure shows the time levels in each run
figTimes = plt.figure(30, facecolor='w')
axTimes = figTimes.add_subplot(1, 1, 1)
plt.xlabel('Year')
axTimesYlabels = []

# =================== plot lag
fig4 = plt.figure(4, facecolor='w')

nrow=3
ncol=1

ax4timelag = fig4.add_subplot(nrow, ncol, 1)
#plt.xlabel('VAF')
plt.ylabel('SLR delay (yr)')
plt.grid()

ax4lagRate = fig4.add_subplot(nrow, ncol, 2, sharex=ax4timelag)
#plt.xlabel('VAF')
plt.ylabel('delay rate (yr (mm SLR)$^{-1}$)')
plt.grid()

#ax4meltDiffRate = fig4.add_subplot(nrow, ncol, 3, sharex=ax4timelag)
#plt.xlabel('VAF')
#plt.ylabel('mean melt rate (m/yr)')
#plt.grid()

ax4VAFrate = fig4.add_subplot(nrow, ncol, 3, sharex=ax4timelag)
plt.xlabel('SLR (mm)')
plt.ylabel('VAF rate (Gt yr$^{-1}$)')
plt.grid()



# =========
print "Done setting up figure axes."
# =========

# --- Define colors for lines ---
colors = []
#colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan']
n150 = sum("amp150" in r for r in runs)
colors.extend( [ cm.autumn(x) for x in np.linspace(0.0, 1.0, n150) ])
n300 = sum("amp300" in r for r in runs)
colors.extend( [ cm.winter(x) for x in np.linspace(0.0, 1.0, n300) ] )
color_index = 0


# ================
# Loop over runs and plot data
# ================
runNumber = 0
for run in runs:
   print "Plotting run: " + run

   thisRun = runData[run]
   # Pull needed data for plotting from this run
   # TODO: replace local variables below in plotting commands with reference to object variables
   yrs=thisRun.yrs
   melt=thisRun.melt
   totalmelt=thisRun.totalmelt
   cumumelt=thisRun.cumumelt
   VAF=thisRun.VAF
   resampYrs=thisRun.resampYrs
   VAFsmoothrate=thisRun.VAFsmoothrate
   GA=thisRun.GA
   GArate=thisRun.GArate
   GLflux = thisRun.GLflux

   # optional rules about coloring
   if "steady" in run:
      lw = 3
      color = 'k'
      #colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, len(runs)) ];  color = colors[color_index];  color_index += 1; lw=1
   elif run=="no-melt":
      lw = 2
      color = [0.7, 0.7, 0.7]  # gray
   else:
      lw = 1
      color = colors[color_index]
      color_index += 1


   # actual plotting begins here
   axMeanMelt.plot(yrs, melt, color=color, linewidth=lw)
   axTotalMelt.plot(yrs, totalmelt, color=color, linewidth=lw)
   axCumuMelt.plot(yrs, cumumelt, color=color, linewidth=lw)

   axVAF.plot(yrs, VAF, label = run, color=color, linewidth=lw)
   #axVAFrate.plot(yrs[1:], VAFrate, color=color, linewidth=lw)
   axVAFrate.plot(resampYrs[1:], VAFsmoothrate, color=color, linewidth=lw)

   axGA.plot(yrs, GA, color=color, linewidth=lw)
   axGArate.plot(resampYrs[1:], GArate, color=color, linewidth=lw)


   # actual plotting begins here for second figure
   ax2MeanMelt.plot(yrs, melt, color=color, linewidth=lw)
   ax2VAF.plot(yrs, VAF, label = run, color=color, linewidth=lw)
   #ax2VAFrate.plot(yrs[1:], VAFrate, color=color, linewidth=lw)
   #ax2VAFrate.plot(resampYrs[1:], VAFsmoothrate, color=color, linewidth=lw)
   ax2VAFrate.plot(resampYrs[windowLength:], VAFsmoothrate[windowLength-1:], color=color, linewidth=lw)
   #ax2VAFrate.plot(resampYrs[:windowLength], VAFsmoothrate[0:windowLength-1], 'x', color=color)  # used to see what values we are throwing away

   ax2GLflux.plot(yrs, GLflux, color=color, linewidth=lw)
   
   ax2floatVol.plot(yrs, thisRun.floatingVol/thisRun.floatingArea, color=color, linewidth=lw)

   # plot run duration
   axTimes.plot(yrs, yrs*0+runNumber, '.')

   if run=='steady':
      # print some stats
      ind1 = np.nonzero(resampYrs<100.0)
      ind2 = np.nonzero(np.logical_and(resampYrs>100.0, resampYrs<200.0))
      print "VAF rate: 1st 100 yrs={}; 2nd 100 yrs={}".format(VAFsmoothrate[ind1].mean(), VAFsmoothrate[ind2].mean())

   # accounting
   axTimesYlabels.append(run)
   runNumber += 1

axVAF.legend(loc='best', ncol=2)

# second ticks for SLR

axSLR=axVAF.twinx()
y1, y2=axVAF.get_ylim()
x1, x2=axVAF.get_xlim()
axSLR.set_ylim(GTtoSL(y1) - GTtoSL(VAF[0]), GTtoSL(y2) - GTtoSL(VAF[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLR.set_ylabel('S.L. equiv. (mm)')
axSLR.set_xlim(x1, x2)

axSLRrate=axVAFrate.twinx()
y1, y2=axVAFrate.get_ylim()
x1, x2=axVAFrate.get_xlim()
axSLRrate.set_ylim(GTtoSL(y1), GTtoSL(y2))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLRrate.set_ylabel('S.L. equiv. (mm/yr)')
axSLRrate.set_xlim(x1, x2)




# Label the progress plot with each run
axTimes.set_yticks(range(runNumber))
axTimes.set_yticklabels(axTimesYlabels) 


# ===================  for figure 2
#handles, labels = ax2VAF.get_legend_handles_labels()
#l1 = ax2VAF.legend(handles[0:1], labels[0:1], loc='lower left', ncol=1)  # control runs
#l2 = ax2VAF.legend(handles[1:], labels[1:], loc='lower center', ncol=2)  # variability runs
#ax2VAF.add_artist(l1)

ax2SLR=ax2VAF.twinx()
y1, y2=ax2VAF.get_ylim()
x1, x2=ax2VAF.get_xlim()
ax2SLR.set_ylim(GTtoSL(y1) - GTtoSL(VAF[0]), GTtoSL(y2) - GTtoSL(VAF[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
ax2SLR.set_ylabel('S.L. equiv. (mm)')
ax2SLR.set_xlim(x1, x2)

ax2SLRrate=ax2VAFrate.twinx()
y1, y2=ax2VAFrate.get_ylim()
x1, x2=ax2VAFrate.get_xlim()
ax2SLRrate.set_ylim(GTtoSL(y1), GTtoSL(y2))
#ax2SLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
ax2SLRrate.set_ylabel('S.L. equiv. (mm/yr)')
ax2SLRrate.set_xlim(x1, x2)

# ================
# Loop over ensemble groups
# ================

# =================== repeat clean version for ensemble plots ===========
# (could make a function to repeat these steps, but they may diverge eventually

fig3 = plt.figure(3, facecolor='w')

nrow=4
ncol=1

# melt forcing
ax3MeanMelt = fig3.add_subplot(nrow, ncol, 1)
plt.xlabel('Year')
plt.ylabel('mean melt (m/yr)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

# VAF
ax3VAF = fig3.add_subplot(nrow, ncol, 2, sharex=ax3MeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF (Gt)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

# VAF diff
#fig9 = plt.figure(9, facecolor='w')
ax3VAFdiff = fig3.add_subplot(nrow, ncol, 4, sharex=ax3MeanMelt)
#ax3VAFdiff = fig9.add_subplot(1, 1, 1)
plt.xlabel('Year', axes=ax3VAFdiff)
#plt.ylabel('VAF difference (%)')
#plt.ylabel('VAF difference (Gt)')
plt.ylabel('SLR difference (mm)', axes=ax3VAFdiff)
#plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()
#plt.xlim((0.0,500.0))

# VAF rate
ax3VAFrate = fig3.add_subplot(nrow, ncol, 3, sharex=ax3MeanMelt)
#plt.xlabel('Year')
plt.ylabel('VAF rate (Gt/yr)')
plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()



# define colors to use
#colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, len(groups)) ]
n150 = sum(1 for g in groups if 'amp150' in g)
colors = [ cm.autumn(x) for x in np.linspace(0.0, 1.0, n150) ]
n300 = sum(1 for g in groups if 'amp300' in g)
colors.extend( [ cm.winter(x) for x in np.linspace(0.0, 1.0, n300) ] )


# ref value
steadyVAF = runData['steady'].VAFsmooth
#steadyVAF = runData['control_run_10pct_D25_kappa8'].VAFsmooth
steadyTimeOnVAFeven = runData['steady'].timeOnVAFeven
steadyMeltOnVAFeven = runData['steady'].meltOnVAFeven
steadyTimeOnGAeven = runData['steady'].timeOnGAeven
steadyMeltOnGAeven = runData['steady'].meltOnGAeven

# add control run
if 'steady' in runData:
   ax3MeanMelt.plot(runData['steady'].resampYrs[windowLength:], runData['steady'].resampMelt[windowLength:], 'k', label='steady')
   ax3VAF.plot(runData['steady'].resampYrs, runData['steady'].VAFsmooth, 'k', label='steady')
   ax3VAFrate.plot(runData['steady'].resampYrs[windowLength:], runData['steady'].VAFsmoothrate[windowLength-1:], 'k', label='steady')

   ax4VAFrate.plot(GTtoSL(runData['steady'].VAFeven) - GTtoSL(runData['steady'].VAF[0]), runData['steady'].VAFrateOnVAFeven, 'k')

groupNumber = 0
for groupName in sorted(groups):  # sorted puts them in alpha order
   print "Plotting group: " + groupName
   group = groups[groupName] # gets the actual dictionary that this group is made of
   nMembers = len(group)

   # Init aggregate stats arrays for this group
   # (could eventually store them in the group if the group was turned into a custom object instead of a dict)
   nEntries = len(group.values()[0].resampYrs)
   yrsGroup = np.zeros((nMembers,))
   meltGroup = np.zeros((nMembers, nEntries))
   VAFgroup = np.zeros((nMembers, nEntries))
   VAFdiffGroup = np.zeros((nMembers, nEntries))
   VAFrateGroup = np.zeros((nMembers, nEntries-1))
   lagGroup = np.zeros((nMembers, len(group.values()[0].VAFeven)))
   meltOnVAFGroup = np.zeros((nMembers, len(group.values()[0].VAFeven)))
   VAFrateOnVAFGroup = np.zeros((nMembers, len(group.values()[0].VAFeven)))
   lagGAGroup = np.zeros((nMembers, len(group.values()[0].GAeven)))
   meltOnGAGroup = np.zeros((nMembers, len(group.values()[0].GAeven)))
#   VAFmean = np.zeros((nMembers,))
#   VAFmin= np.zeros((nMembers,))
#   VAFmax= np.zeros((nMembers,))

   # Loop through group members
   runNumber = 0
   for run in group:

      thisRun = group[run]

      # use the resampled versions so they all are on the same axis
      yrsGroup = thisRun.resampYrs  # (only need this once)
      meltGroup[runNumber, :] = thisRun.resampMelt
      VAFgroup[runNumber, :] = thisRun.VAFsmooth
      #VAFdiffGroup[runNumber, :] = thisRun.VAFsmooth - steadyVAF # this version plots difference from control
      VAFdiffGroup[runNumber, :] = GTtoSL(thisRun.VAFsmooth-steadyVAF)  # this version plots difference from control as SLR
    #  VAFdiffGroup[runNumber, :] = (thisRun.VAFsmooth - steadyVAF)/(steadyVAF-steadyVAF[0])*100.0 # this version plots difference from control as %
      VAFrateGroup[runNumber, :] = thisRun.VAFsmoothrate
      meltOnVAFGroup[runNumber, :] = thisRun.meltOnVAFeven 
      meltOnGAGroup[runNumber, :] = thisRun.meltOnGAeven 

      #VAFevenGroup = thisRun.VAFeven # only need this once)
      VAFevenGroup = GTtoSL(thisRun.VAFeven-thisRun.VAF[0]) # only need this once - as SLR
      lagGroup[runNumber,:] = thisRun.timeOnVAFeven - steadyTimeOnVAFeven
      VAFrateOnVAFGroup[runNumber,:] = thisRun.VAFrateOnVAFeven
#      ax4timelag.plot(thisRun.VAFeven, lagGroup[runNumber,:], ':', color=colors[groupNumber], linewidth=0.5)
#      ax4VAFrate.plot(thisRun.VAFeven, thisRun.VAFrateOnVAFeven, color=colors[groupNumber], linewidth=0.4)
      #ax4meltDiffRate.plot(thisRun.VAFeven, thisRun.meltOnVAFeven - steadyMeltOnVAFeven, color=colors[groupNumber])
      #ax4meltDiffRate.plot(thisRun.VAFeven, thisRun.meltOnVAFeven , color=colors[groupNumber], linewidth=0.4)

      GAevenGroup = thisRun.GAeven # only need this once)
      lagGAGroup[runNumber,:] = thisRun.timeOnGAeven - steadyTimeOnGAeven
      #ax4timelag.plot(thisRun.GAeven, lagGAGroup[runNumber,:], ':', color=colors[groupNumber], linewidth=0.5)
      #ax4VAFrate.plot(thisRun.GAeven, thisRun.VAFrateOnVAFeven, color=colors[groupNumber])

      runNumber += 1
   
   # melt plot
   ax3MeanMelt.plot(yrsGroup[windowLength:], meltGroup.mean(0)[windowLength:], '-', color = colors[groupNumber], label=groupName)
   ax3MeanMelt.plot(yrsGroup[windowLength:], meltGroup.max(0)[windowLength:], '--', color = colors[groupNumber], linewidth=0.5)
   ax3MeanMelt.plot(yrsGroup[windowLength:], meltGroup.min(0)[windowLength:], '--', color = colors[groupNumber], linewidth=0.5)

   # VAF plot
   ax3VAF.plot(yrsGroup, VAFgroup.mean(0), '-', color = colors[groupNumber], label=groupName)
   ax3VAF.plot(yrsGroup, VAFgroup.max(0), '--', color = colors[groupNumber], linewidth=0.5)
   ax3VAF.plot(yrsGroup, VAFgroup.min(0), '--', color = colors[groupNumber], linewidth=0.5)

   # VAF diff plot
   ax3VAFdiff.plot(yrsGroup, VAFdiffGroup.mean(0), '-', color = colors[groupNumber], label=groupName)
   ax3VAFdiff.plot(yrsGroup, VAFdiffGroup.max(0), '--', color = colors[groupNumber], linewidth=0.5)
   ax3VAFdiff.plot(yrsGroup, VAFdiffGroup.min(0), '--', color = colors[groupNumber], linewidth=0.5)

   # VAF rate plot
   ax3VAFrate.plot(yrsGroup[windowLength:], VAFrateGroup.mean(0)[windowLength-1:], '-', color = colors[groupNumber], label=groupName)
   ax3VAFrate.plot(yrsGroup[windowLength:], VAFrateGroup.max(0)[windowLength-1:], '--', color = colors[groupNumber], linewidth=0.5)
   ax3VAFrate.plot(yrsGroup[windowLength:], VAFrateGroup.min(0)[windowLength-1:], '--', color = colors[groupNumber], linewidth=0.5)

   # plot DELAY
   ax4timelag.plot(VAFevenGroup, lagGroup.mean(0), '-', color = colors[groupNumber], label=groupName)
   ax4timelag.plot(VAFevenGroup, lagGroup.max(0), '--', color = colors[groupNumber], linewidth=0.5)
   ax4timelag.plot(VAFevenGroup, lagGroup.min(0), '--', color = colors[groupNumber], linewidth=0.5)
   #ax4timelag.plot(VAFevenGroup, lagGroup.max(0) - lagGroup.min(0), '-', color = colors[groupNumber], linewidth=1)
   for t in np.arange(0.0, 500.0, 20.0):
      idx =  np.abs(steadyTimeOnVAFeven - t).argmin()
      ax4timelag.plot(VAFevenGroup[idx], lagGroup.mean(0)[idx], 'o', color = colors[groupNumber], markersize=3)
   for t in np.arange(0.0, 500.0, 100.0):
      idx =  np.abs(steadyTimeOnVAFeven - t).argmin()
      ax4timelag.plot(VAFevenGroup[idx], lagGroup.mean(0)[idx], 'o', color = colors[groupNumber], markersize=5)

   # plot lag rate
   step=120
   ax4lagRate.plot(VAFevenGroup[step/2:-1*step/2], (lagGroup.mean(0)[:-1*step]-lagGroup.mean(0)[step:]) / (VAFevenGroup[:-1*step]-VAFevenGroup[step:]), '-', color=colors[groupNumber])

#   # plot melt on VAF to go with lag
##   ax4meltDiffRate.plot(VAFevenGroup, meltOnVAFGroup.mean(0)- steadyMeltOnVAFeven, color=colors[groupNumber])
#   ax4meltDiffRate.plot(VAFevenGroup, meltOnVAFGroup.mean(0), color=colors[groupNumber])
#   ax4meltDiffRate.plot(VAFevenGroup, meltOnVAFGroup.max(0), '-', color=colors[groupNumber], linewidth=0.5)
#   ax4meltDiffRate.plot(VAFevenGroup, meltOnVAFGroup.min(0), '-', color=colors[groupNumber], linewidth=0.5)
##   ax4meltDiffRate.plot(VAFevenGroup, np.convolve(np.ones(511,'d')/511.0,meltOnVAFGroup.mean(0),mode='same'), color=colors[groupNumber], linewidth=2)  # this line plots a boxcar filtered version
#   ax4meltDiffRate.plot(VAFevenGroup, steadyMeltOnVAFeven, color='k')
   print groupName, "variability mean melt rate: ", meltOnVAFGroup.mean(0).mean(), " control mean melt rate: ", steadyMeltOnVAFeven.mean(), " ratio=", meltOnVAFGroup.mean(0).mean()/steadyMeltOnVAFeven.mean()

   # plot VAF rate on VAF to go with lag
   ax4VAFrate.plot(VAFevenGroup, VAFrateOnVAFGroup.mean(0), '-', color = colors[groupNumber], label=groupName)
   ax4VAFrate.plot(VAFevenGroup, VAFrateOnVAFGroup.max(0), '-', color = colors[groupNumber], linewidth=0.5)
   ax4VAFrate.plot(VAFevenGroup, VAFrateOnVAFGroup.min(0), '-', color = colors[groupNumber], linewidth=0.5)
#   # plot lag based on GA
#   ax4timelag.plot(GAevenGroup, lagGAGroup.mean(0), '-', color = colors[groupNumber], label=groupName)
#   ax4timelag.plot(GAevenGroup, lagGAGroup.max(0), '--', color = colors[groupNumber], linewidth=0.5)
#   ax4timelag.plot(GAevenGroup, lagGAGroup.min(0), '--', color = colors[groupNumber], linewidth=0.5)
#   #ax4timelag.plot(VAFevenGroup, lagGroup.max(0) - lagGroup.min(0), '-', color = colors[groupNumber], linewidth=1)
#   for t in np.arange(0.0, 500.0, 20.0):
#      idx =  np.abs(steadyTimeOnGAeven - t).argmin()
#      ax4timelag.plot(GAevenGroup[idx], lagGAGroup.mean(0)[idx], 'o', color = colors[groupNumber], markersize=3)
#   for t in np.arange(0.0, 500.0, 100.0):
#      idx =  np.abs(steadyTimeOnGAeven - t).argmin()
#      ax4timelag.plot(GAevenGroup[idx], lagGAGroup.mean(0)[idx], 'o', color = colors[groupNumber], markersize=5)
#
#   # plot lag rate
#   step=30
#   ax4lagRate.plot(GAevenGroup[step/2:-1*step/2], (lagGAGroup.mean(0)[:-1*step]-lagGAGroup.mean(0)[step:]) / (GAevenGroup[:-1*step]-GAevenGroup[step:]), '-', color=colors[groupNumber])
#
##   ax4meltDiffRate.plot(VAFevenGroup, meltOnVAFGroup.mean(0)- steadyMeltOnVAFeven, color=colors[groupNumber])
#   ax4meltDiffRate.plot(GAevenGroup, meltOnGAGroup.mean(0), color=colors[groupNumber])
#   ax4meltDiffRate.plot(GAevenGroup, steadyMeltOnGAeven, color='k')


   groupNumber += 1

# show legend   
#legend = ax3VAF.legend(loc='lower left')
# or as multiple columns
handles, labels = ax3VAF.get_legend_handles_labels()
l1 = ax3VAF.legend(handles[0:1], labels[0:1], loc='lower left', ncol=1)  # control runs
l2 = ax3VAF.legend(handles[1:], labels[1:], loc='lower center', ncol=2)  # variability runs
ax3VAF.add_artist(l1)


axSLR=ax3VAF.twinx()
y1, y2=ax3VAF.get_ylim()
x1, x2=ax3VAF.get_xlim()
axSLR.set_ylim(GTtoSL(y1) - GTtoSL(runData['steady'].VAFsmooth[0]), GTtoSL(y2) - GTtoSL(runData['steady'].VAFsmooth[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLR.set_ylabel('S.L. equiv. (mm)')
axSLR.set_xlim(x1, x2)

#axSLR=ax3VAFdiff.twinx()
#y1, y2=ax3VAFdiff.get_ylim()
#x1, x2=ax3VAFdiff.get_xlim()
#axSLR.set_ylim(GTtoSL(y1) , GTtoSL(y2) )
##axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
#axSLR.set_ylabel('S.L. equiv. (mm)')
#axSLR.set_xlim(x1, x2)



plt.show()
