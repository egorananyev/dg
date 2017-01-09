####!/usr/bin/arch -i386 /usr/bin/python # -*- coding: utf-8 -*-
"""
Drifting Gratings
2016-07-29
"""

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, event, gui
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np
import pandas as pd
from datetime import datetime
import os, shutil, itertools  # handy system and path functions

# Initiating the keyboard
from psychopy.iohub import launchHubServer
io = launchHubServer()
kb_device = io.devices.keyboard

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))

# ====================================================================================
## Initial variables.
et = 1
expName = 'dg'
# Window circles (specified in degrees of visual angles [dva]):
#winSz = 7.2 # 5.03; calculated as 5/x=sqrt(2)/2 => x=10/sqrt(2)
winOffX = 4.25 # 6 # 5.62
winOffY = 3.5 # 5.5 (3.5cm ~= 124px)
winThickness = 2 # in pixels
fdbkLen = .25 # the length of the feedback line, in degrees
fdbkThick = 4 # the tickness of the feedback line, in pixels
ISIduration = 1
grtSize = 256 # size of 256 is 71mm, or 7.2dova
# color mask:
colMaskAlpha = .2 # default alpha value
colMaskMagentaSat = .25 # default magenta saturation 
# Dimensions:
###### 7.2dova = 71mm = 256px; 475x296mm, 563mm viewing dist ######
dr = (1680,1050) # display resolution in px
dd = (29.5,16.6) # display dimensions in cm
ds = 50+2.5+3.5 #49.5 # distance to screen in cm
trialNfb = False # do we give the trial number feedback?
nFrames = 60 # number of frames per sequence

# ====================================================================================
# Converter functions:
def cm2px(cm,dr=dr,dd=dd):
    px = int(cm*(dr[0]/dd[0]))
    return px
def px2cm(px,dr=dr,dd=dd):
    cm = px/(dr[0]/dd[0])
    return cm
def cm2dg(cm,ds=ds):
    dg = np.degrees(np.arctan(cm/ds))
    return dg
def dg2cm(dg,ds=ds):
    cm = ds*np.tan(np.radians(dg))
    return cm
def px2dg(px,cm2dg=cm2dg,px2cm=px2cm):
    dg = cm2dg(px2cm(px))
    return dg
def dg2px(dg,cm2px=cm2px,dg2cm=dg2cm):
    px = int(cm2px(dg2cm(dg)))
    return px

# ====================================================================================
# Converting win dimensions to pixels
#winSz = dg2px(winSz)
winSz = grtSize + 2
winOffX = dg2px(winOffX)
winOffY = dg2px(winOffY)
fdbkLen = dg2px(fdbkLen)
posCentL = [-winOffX, winOffY]
posCentR = [winOffX, winOffY]
print winSz 
print posCentL 
print posCentR 
print fdbkLen 

# ====================================================================================
# Eye tracking initialization

if et:
    import pylink as pl
    #cp = (0.4,0.4) # calibration proportion
    cd = 32

    eyeLink = ("100.1.1.1")

    displayInfo = pl.getDisplayInformation()
    print displayInfo.width, displayInfo.height
    screenCenter = (int(dr[0]/2), int(dr[1]/2))
    calScreenCenter = (int(screenCenter[0]+winOffX),
                    int(screenCenter[1]-winOffY))
    calTargDist = int(winSz/3)
    calTarg1 = calScreenCenter
    calTarg2 = (int(calScreenCenter[0]-calTargDist), int(calScreenCenter[1]))
    calTarg3 = (int(calScreenCenter[0]+calTargDist), int(calScreenCenter[1]))

    def elEndRec(el):
        # Ends the recording; adds 100ms to catch final events
        pl.endRealTimeMode()
        pl.pumpDelay(100)
        el.stopRecording()

    def eyeTrkInit (dr):
        el = pl.EyeLink()
        # sending the screen dimensions to the eye tracker:
        el.sendCommand('screen_pixel_coords = 0 0 %d %d' %dr)
        el.sendMessage('DISPLAY_COORDS 0 0 %d %d' %dr)
        el.sendCommand('generate_default_targets = NO')
        el.sendCommand('calibration_targets = %d,%d %d,%d %d,%d' % (
                       calTarg1[0], calTarg1[1],
                       calTarg2[0], calTarg2[1],
                       calTarg3[0], calTarg3[1]) )
        el.sendCommand('validation_targets = %d,%d %d,%d %d,%d' % (
                       calTarg1[0], calTarg1[1],
                       calTarg2[0], calTarg2[1],
                       calTarg3[0], calTarg3[1]) )
        # parser configuration 1 corresponds to high sensitivity to saccades:
        el.sendCommand('select_parser_configuration 1')
        # turns off "scenelink camera stuff", i.e., doesn't record the ET video
        el.sendCommand('scene_camera_gazemap = NO')
        # converting pupil area to diameter
        el.sendCommand('pupil_size_diameter = %s'%('YES'))
        return(el)
    el = eyeTrkInit(dr)
    print 'Finished initializing the eye tracker.'

    def eyeTrkCalib (el=el,dr=dr,cd=cd):
        # "opens the graphics if the display mode is not set"
        pl.openGraphics(dr,cd)
        pl.setCalibrationColors((255,255,255),(0,177,177))
        pl.setTargetSize(10, 5) 
        pl.setCalibrationSounds("","","")
        el.setCalibrationType('H3')
        pl.setDriftCorrectSounds("","off","off")
        el.disableAutoCalibration()
        el.doTrackerSetup()
        el.drawCalTarget(calTarg1)
        el.drawCalTarget(calTarg2)
        el.drawCalTarget(calTarg3)
        pl.closeGraphics()
        el.setOfflineMode()

# ====================================================================================
# Store info about the experiment session and run:
expInfo = {u'session': u'', u'run': u'', u'participant': u'', u'condition': u'rn-test'}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName) # dialogue box
if dlg.OK == False: core.quit()  # user pressed cancel
timeNow = datetime.now()
expInfo['time'] = datetime.now().strftime('%Y-%m-%d_%H%M')
expInfo['expName'] = expName
curCond = expInfo['condition'] # condition: report/no report (r/n) attnTask/no attnTask (a/n)

# Setup the window
win = visual.Window(size=dr, fullscr=True, screen=0, allowGUI=False, 
      allowStencil=False, color='grey', blendMode='avg', useFBO=True, units='pix')
# store frame rate of monitor if we can measure it successfully:
frameRate=win.getActualFrameRate()
if frameRate!=None:
    frameDur = 1.0/round(frameRate)
else:
    frameDur = 1.0/60.0 # couldn't get a reliable measure so guess

# ====================================================================================
# Eye-tracking setup

if et:
    def etSetup(el=el,dr=dr,cd=cd):
        blockLabel=visual.TextStim(win, text="Press spacebar", pos=posCentR,
                                   color="white", bold=True, alignHoriz="center",
                                   height=0.5)
        notdone=True
        while notdone:
            blockLabel.draw()
            win.flip()
            keySpace = event.getKeys(keyList=['escape','space'])
            if 'space' in keySpace:
                print 'spacebar pressed'
                eyeTrkCalib()
                win.winHandle.activate()
                print '///Finished calibration///'
                notdone=False
            elif 'escape' in keySpace:
                print 'procedure terminated'
                notdone=False
    etSetup()

    def drCor(el=el,dr=dr,cd=cd):
        pl.openGraphics(dr,cd)
        el.doDriftCorrect(calScreenCenter[0], calScreenCenter[1], 1, 0)
        pl.closeGraphics()
        print '///Finished drift correction///'

# ====================================================================================

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
dataDir = '..' + os.sep + 'data'
fileName = '%s_p%s_c%s_s%s_r%s_%s' %(expName, expInfo['participant'], curCond, 
    expInfo['session'], expInfo['run'],expInfo['time'])
filePath = dataDir + os.sep + fileName
print filePath

if et:
    edfFileName = 'data.edf'
    el.openDataFile(edfFileName)
    el.sendCommand("file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,\
                    MESSAGE,BUTTON,INPUT")
    el.sendCommand("file_sample_data  = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS,\
                    HTARGET,INPUT")
    print '///set up the EDF file for eye-tracking///'

# Condition-related variables
conditionsFilePath = 'cond-files'+os.sep+'cond-'+expName+'-'+curCond+'.csv'
print conditionsFilePath
os.chdir(_thisDir)

# ====================================================================================

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Initialize components for Routine "instructions"
instructionsClock = core.Clock()
instrTextNN = 'Please passively observe the stimuli; no need to press anything.'
instrTextRN = 'Please report the direction of the currently visible drifting stimulus when it is fully visible.'
instrTextNA = 'Please report the pitch of the tone by pressing arrow keys ("up" for high, "down" for low).'
if curCond == 'nn':
    curInstrText = instrTextNN
elif curCond == 'na':
    curInstrText = instrTextNA
elif curCond == 'rn' or curCond == 'rn-test':
    curInstrText = instrTextRN
else:
    print('ERROR: condition not recognized.')
instrTextL = visual.TextStim(win, text=curInstrText, font='Cambria',
                             pos=posCentL, height=dg2px(.45), wrapWidth=dg2px(4.5),
                             color='white', alignHoriz='center')
instrTextR = visual.TextStim(win, text=curInstrText, font='Cambria',
                             pos=posCentR, height=dg2px(.45), wrapWidth=dg2px(4.5),
                             color='white', alignHoriz='center')

# Initialize components for Routine "trial"
trialClock = core.Clock()
moveClock = core.Clock()
maskMoveClock = core.Clock()
ISI = core.StaticPeriod(win=win, screenHz=frameRate, name='ISI')
# circular windows:
winL = visual.Polygon(win, edges=36, size=[winSz, winSz], pos=posCentL,
                      lineWidth=winThickness, lineColor='white')
winR = visual.Polygon(win, edges=36, size=[winSz, winSz], pos=posCentR,
                      lineWidth=winThickness, lineColor='white')
# gabor patches:
gabL = visual.GratingStim(win, tex='sin', mask='circle', size=[winSz, winSz],
                          pos=posCentL)
gabR = visual.GratingStim(win, tex='sin', mask='circle', size=[winSz, winSz],
                          pos=posCentR)
# color masks:
#colMaskL = visual.GratingStim(win, size=[grtSize, grtSize], pos=posCentL, opacity=colMaskAlpha,
#                              colorSpace='hsv', mask='circle')
colMaskL = visual.Polygon(win, edges=36, size=[grtSize, grtSize], pos=posCentL, opacity=colMaskAlpha)
colMaskR = visual.Polygon(win, edges=36, size=[grtSize, grtSize], pos=posCentR, opacity=colMaskAlpha)
# direction feedback:
dirFdbkL = visual.Line(win, start=[0,0], end=[0,0], lineColor='white',
                        lineWidth=fdbkThick)
dirFdbkR = visual.Line(win, start=[0,0], end=[0,0], lineColor='white',
                        lineWidth=fdbkThick)
# Trial number feedback:
if trialNfb:
    trialNfbText = visual.TextStim(win=win, text='', font='Cambria', 
                                   pos=(0,0), height=dg2px(.55), wrapWidth=dg2px(4.5),
                                   color='white')
# pause text:
pauseTextL = visual.TextStim(win, text='Press Spacebar to continue', font='Cambria',
                             alignHoriz='center', pos=posCentL, height=dg2px(.7),
                             wrapWidth=dg2px(3), color='white')
pauseTextR = visual.TextStim(win, text='Press Spacebar to continue', font='Cambria',
                             alignHoriz='center', pos=posCentR, height=dg2px(.7),
                             wrapWidth=dg2px(3), color='white')

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

#------Prepare to start Routine "instructions"-------
t = 0
instructionsClock.reset()  # clock 
frameN = -1
# update component parameters for each repeat
instrKey = event.BuilderKeyResponse()  # create an object of type KeyResponse
instrKey.status = NOT_STARTED
# keep track of which components have finished
instructionsComponents = []
instructionsComponents.append(instrTextL)
instructionsComponents.append(instrTextR)
instructionsComponents.append(instrKey)
for thisComponent in instructionsComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED


# ====================================================================================
# Setting up the conditions:
condList = data.importConditions(conditionsFilePath)
conds = []
commonNTrials = []
for thisCondition in condList:
    nTrials = thisCondition['trialN']
    # print 'Number of trials in this condition: ' + str(nTrials)
    conds.append(thisCondition)
    commonNTrials = nTrials

# An empty data set for storing behavioural responses:
behResp = []
    
# Printing the attributes of the conds:  
print commonNTrials
trials = data.TrialHandler(conds, commonNTrials, extraInfo=expInfo)

# Creating a copy of the Conditions file for book-keeping and analyses:
if not os.path.exists(filePath):
    os.makedirs(filePath)
shutil.copyfile(conditionsFilePath, filePath + os.sep + 
                os.path.basename(conditionsFilePath))
dataFileName = filePath + os.sep + fileName + '.csv'

# ====================================================================================
# Various functions for use in trials:

def drawFdbkAngle(dirFdbk, lr, angle, winOffX=winOffX, winOffY=winOffY, winSz=winSz):
    lineStart = [int( lr*winOffX + np.cos(np.radians(angle))*winSz/2 ),
                 int( winOffY + np.sin(np.radians(angle))*winSz/2 )]
    lineEnd = [int( lr*winOffX + np.cos(np.radians(angle)) * (winSz/2 + fdbkLen) ),
               int( winOffY + np.sin(np.radians(angle)) * (winSz/2 + fdbkLen) )]
    dirFdbk.start = lineStart
    dirFdbk.end = lineEnd
    return dirFdbk

# ====================================================================================

#-------Start Routine "instructions"-------
continueRoutine = True
while continueRoutine:
    # get current time
    t = instructionsClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *instrText* updates
    if t >= 0.0 and instrTextL.status == NOT_STARTED:
        # keep track of start time/frame for later
        instrTextL.tStart = t  # underestimates by a little under one frame
        instrTextL.frameNStart = frameN  # exact frame index
        instrTextL.setAutoDraw(True)
        instrTextR.tStart = t  # underestimates by a little under one frame
        instrTextR.frameNStart = frameN  # exact frame index
        instrTextR.setAutoDraw(True)
    
    # *instrKey* updates
    if t >= 0.0 and instrKey.status == NOT_STARTED:
        # keep track of start time/frame for later
        instrKey.tStart = t  # underestimates by a little under one frame
        instrKey.frameNStart = frameN  # exact frame index
        instrKey.status = STARTED
        # keyboard checking is just starting
        event.clearEvents(eventType='keyboard')
        winL.setAutoDraw(True)
        winR.setAutoDraw(True)
    if instrKey.status == STARTED:
        theseKeys = event.getKeys()
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            # a response ends the routine
            continueRoutine = False
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineTimer.reset()  # if we abort early the non-slip timer needs reset
        break
    continueRoutine = False  # reverts to True if at least 1 component still running
    for thisComponent in instructionsComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()
    else:  # this Routine was not non-slip safe so reset non-slip timer
        routineTimer.reset()

#-------Ending Routine "instructions"-------
for thisComponent in instructionsComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

#win.winHandle.minimize()
#win.flip()
#drCor()
#win.winHandle.maximize()
#win.flip()
#win.winHandle.activate()

# ====================================================================================
# Initiating the trial loop

nDone=0
for thisTrial in trials:
    print '===new=trial==='
    nDone += 1
    if trialNfb:
        trialNfbText.text = str(nDone) + '/' + str(trials.nTotal)
    print 'trial#' + str(nDone) + '/' + str(trials.nTotal)

    ## Setting up trial variables

    # Setting static Gabor parameters here:
    szL = thisTrial['szL']
    szR = thisTrial['szR']
    sfL = thisTrial['sfL']
    sfR = thisTrial['sfR']
    gabL.sf = sfL
    gabR.sf = sfR
    gabL.size = szL
    gabR.size = szR

    # Motion parameters for the trial:
    dirL = thisTrial['dirL']
    dirR = thisTrial['dirR']
    vL = thisTrial['vL']
    vR = thisTrial['vR']

    # Other variables:
    colL = thisTrial['colL']
    colR = thisTrial['colR']
    trialT = thisTrial['trialT'] # -win.monitorFramePeriod*0.75

    # Print out trial info based on the nature of the experiment:
    print 'dirL=' + str(dirL) + '; dirR=' + str(dirR)
    print 'vL=' + str(vL) + '; vR=' + str(vR)

    # Color, if any:
    #colorEither = [[150,1,1],[330,colMaskMagentaSat,1]] # green and magenta
    colorEither = [[0,1,0],[1,0,0]] # green and red
    if thisTrial['colL'] == 'rand': # picking one at random
        colorPick = np.random.permutation(colorEither)
        colorL = colorPick[0]
        colorR = colorPick[1]
    elif thisTrial['colL'] == thisTrial['colR']:
        colorL = [0,0,0] # greyscale (black)
        colorR = [0,0,0]
    else:
        if thisTrial['colL'] == 'green':
            colorL = colorEither[0] # green
            colorR = colorEither[1] # red
        else:
            colorL = colorEither[1] # red
            colorR = colorEither[0] # green
    print 'colorL = ' + str(colorL) + '; colorR = ' + str(colorR)
    
    # Creating an empty matrix for keeping the behavioural responses:
    # enable continuous tracking:
    behRespTrial = np.empty([1, trialT*nFrames]) 
    behRespTrial[:] = np.NAN
        
    # Using the mask to assign both the greyscale values and the mask for our color masks:
    if not colorL == colorR: # same colors mean no color mask
        #colMaskL.tex = (curMaskL + 1) / 2
        colMaskL.fillColor = colorL
        #colMaskL.mask = curMaskL
        #colMaskR.tex = (curMaskR + 1) / 2
        colMaskR.fillColor = colorR
        #colMaskR.mask = curMaskR

    #------Prepare to start Routine "trial"-------
    t = 0
    trialClock.reset()  # clock 
    frameN = -1
    tMaskMove = 0
    elStopped = False
    key_pause = False
    behRespRecorded = False
    someKeyPressed = False # to prevent recording key releases at trial beginning
    # update component parameters for each repeat
    key_arrow = event.BuilderKeyResponse()  # create an object of type KeyResponse
    key_arrow.status = NOT_STARTED
    # keep track of which components have finished
    trialComponents = []
    trialComponents.append(winL)
    trialComponents.append(winR)
    trialComponents.append(dirFdbkL)
    trialComponents.append(dirFdbkR)
    trialComponents.append(key_arrow)
    trialComponents.append(pauseTextL)
    trialComponents.append(pauseTextR)
    trialComponents.append(ISI)
    for thisComponent in trialComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # ////////////////////////////////////////////////////////////////////////////////
    #win.winHandle.minimize()
    #drCor(el,dr,cd)
    #win.winHandle.maximize()
    #win.winHandle.activate()
    if et:
        el.sendMessage("TRIALID " + str(nDone))
        trialStartStr = datetime.now().strftime('%Y-%m-%d_%H%M%S')
        el.sendMessage("TIMESTAMP " + trialStartStr)
        el.setOfflineMode()
        pl.msecDelay(50) 
        error = el.startRecording(1,1,1,1)
    # ////////////////////////////////////////////////////////////////////////////////
    
    #-------Start Routine "trial"-------
    continueRoutine = True
    while continueRoutine:
        # get current time
        t = trialClock.getTime()
        frameN = frameN + 1 # number of completed frames (0 is the first frame)
        # update/draw components on each frame
        
        # *winL* updates
        if winL.status == NOT_STARTED:
            # keep track of start time/frame for later
            winL.tStart = t  # underestimates by a little under one frame
            winL.frameNStart = frameN  # exact frame index
            winL.setAutoDraw(True)
            winL.status = STARTED
        
        # *winR* updates
        if winR.status == NOT_STARTED:
            # keep track of start time/frame for later
            winR.tStart = t  # underestimates by a little under one frame
            winR.frameNStart = frameN  # exact frame index
            winR.setAutoDraw(True)
            winR.status = STARTED

        if trialNfb:
            trialNfbText.draw()

        # stimulus presentation:
        if t < trialT:
            # grating motion:
            gabL.phase = gabL.phase + np.cos(dirL) * vL / 60
            gabR.phase = gabR.phase + np.cos(dirR) * vR / 60
            # drawing the grating and the mask:
            gabL.draw()
            gabR.draw()
            colMaskL.draw()
            colMaskR.draw()
        
        # *key_arrow* updates
        if key_arrow.status == NOT_STARTED:
            # keep track of start time/frame for later
            key_arrow.tStart = t  # underestimates by a little under one frame
            key_arrow.frameNStart = frameN  # exact frame index
            key_arrow.status = STARTED
            # keyboard checking is just starting
            key_arrow.clock.reset()  # now t=0
            event.clearEvents(eventType='keyboard')
            kb_device.clearEvents()
        # registering the response continuously:
        if key_arrow.status == STARTED and t < trialT:
            thesePresses = kb_device.getPresses(keys=['left','right','up','down'])
            theseReleases = kb_device.getReleases(keys=['left','right','up','down'])
            if len(thesePresses) > 0:
                dirFdbkL.setAutoDraw(True)
                dirFdbkR.setAutoDraw(True)
                keyPressFN = frameN
                someKeyPressed = True
                if 'left' in thesePresses:
                    print '"left" key is pressed'
                    drawFdbkAngle(dirFdbkL, -1, 180)
                    drawFdbkAngle(dirFdbkR, 1, 180)
                    whichKeyPressed = 'left' # only needed for final key press
                elif 'right' in thesePresses:
                    print '"right" key is pressed'
                    drawFdbkAngle(dirFdbkL, -1, 0)
                    drawFdbkAngle(dirFdbkR, 1, 0)
                    whichKeyPressed = 'right'
                elif 'up' in thesePresses:
                    print '"up" key is pressed'
                    drawFdbkAngle(dirFdbkL, -1, 90)
                    drawFdbkAngle(dirFdbkR, 1, 90)
                    whichKeyPressed = 'up'
                elif 'down' in thesePresses:
                    print '"down" key is pressed'
                    drawFdbkAngle(dirFdbkL, -1, 270)
                    drawFdbkAngle(dirFdbkR, 1, 270)
                    whichKeyPressed = 'down'
            if len(theseReleases) > 0 and someKeyPressed:
                dirFdbkL.setAutoDraw(False)
                dirFdbkR.setAutoDraw(False)
                print '...released after ' + \
                      str(np.around((frameN-keyPressFN)/60,2)) + 's'
                someKeyPressed = False
                if 'left' in theseReleases:
                    behRespTrial[0,keyPressFN:frameN+1] = 180 # left
                elif 'right' in theseReleases:
                    behRespTrial[0,keyPressFN:frameN+1] = 0 # right
                elif 'up' in theseReleases:
                    behRespTrial[0,keyPressFN:frameN+1] = 90 # right
                elif 'down' in theseReleases:
                    behRespTrial[0,keyPressFN:frameN+1] = 270 # down

        if t > trialT and not elStopped:
            # stopping eye-tracking recording:
            if et:
                elEndRec(el)
                elStopped = True

        # pause text and data exporting
        if not key_pause and t>trialT:
            if not behRespRecorded: # a flag for data recording
                # Make sure to record the release of a key at trial end
                if someKeyPressed:
                    if whichKeyPressed == 'left':
                        behRespTrial[0,keyPressFN:(trialT*nFrames)] = 180
                    if whichKeyPressed == 'right':
                        behRespTrial[0,keyPressFN:(trialT*nFrames)] = 0
                    if whichKeyPressed == 'up':
                        behRespTrial[0,keyPressFN:(trialT*nFrames)] = 90
                    if whichKeyPressed == 'down':
                        behRespTrial[0,keyPressFN:(trialT*nFrames)] = 270
                    dirFdbkL.setAutoDraw(False)
                    dirFdbkR.setAutoDraw(False)
                    print '...keypress forced end after ' + \
                        str(np.around(((trialT*nFrames)-keyPressFN)/60,2)) + 's'
                    print 'recorded post-trial response'
                # Recording the responses:
                pauseTextL.setAutoDraw(True)
                pauseTextR.setAutoDraw(True)
                behRespRecorded = True
            if 'space' in event.getKeys(keyList=['space']):
                # Computing and recording predominance:
                nNa = np.count_nonzero(np.isnan(behRespTrial))
                nf000 = np.count_nonzero(behRespTrial==0)
                nf090 = np.count_nonzero(behRespTrial==90)
                nf180 = np.count_nonzero(behRespTrial==180)
                nf270 = np.count_nonzero(behRespTrial==270)
                dfContT = behRespTrial # continuous response for this trial
                dT = pd.DataFrame({'expName': expName,
                                'time': expInfo['time'], 'participant': expInfo['participant'],
                                'session': expInfo['session'], 'run': expInfo['run'],
                                'trialN': nDone, 'dirL': dirL, 'dirR': dirR,
                                'vL': vL, 'vR': vR, 'szL': szL, 'szR': szR,
                                'sfL': sfL, 'sfR': sfR, 'colL': 
                                thisTrial['colL'], 'colR': thisTrial['colR'], 
                                'sat': str(colMaskMagentaSat), 'trialT': trialT, 'nFrames': nFrames, 
                                'nNa': nNa,
                                'nf000': nf000, 'nf090': nf090, 'nf180': nf180, 'nf270': nf270,
                                'pd000': [nf000 / (trialT * nFrames)],
                                'pd090': [nf090 / (trialT * nFrames)],
                                'pd180': [nf180 / (trialT * nFrames)],
                                'pd270': [nf270 / (trialT * nFrames)] })
                # to preserve the column order:
                dataCols = ['expName', 'time', 'participant', 'session', 'run', 'trialN',
                            'dirL', 'dirR', 'vL', 'vR', 'szL', 'szR', 'sfL', 'sfR', 
                            'colL', 'colR', 'sat', 'trialT', 'nFrames', 'nNa',
                            'nf000', 'nf090', 'nf180', 'nf270', 'pd000', 'pd090', 'pd180', 'pd270']
                if nDone == 1:
                    df = dT
                    dfCont = dfContT
                else:
                    df = pd.concat([df,dT])
                    #dfCont = pd.concat([dfCont,dfContT])
                    dfCont = np.concatenate((dfCont,dfContT),axis=0)
                # Recording the data to a csv file:
                df.to_csv(dataFileName, index=False, columns=dataCols)
                #dfCont.to_csv(filePath + os.sep + fileName + '_cont.csv', index=False)
                np.savetxt(filePath + os.sep + fileName + '_cont.csv', dfCont, delimiter=',')
                print 'wrote the data set to ' + dataFileName
                print 'spacebar pressed - continuing to the next trial'
                pauseTextL.setAutoDraw(False)
                pauseTextR.setAutoDraw(False)
                key_pause = True

        # *ISI* period
        if ISI.status == NOT_STARTED and t>=trialT and key_pause:
            # keep track of start time/frame for later
            ISI.tStart = t  # underestimates by a little under one frame
            ISI.frameNStart = frameN  # exact frame index
            ISI.start(ISIduration)
        #one frame should pass before updating params and completing
        elif ISI.status == STARTED and t >= (ISI.tStart + ISIduration): 
            ISI.complete() #finish the static period
            continueRoutine = False
        
        # check if all components have finished
        # a component has requested a forced-end of Routine:
        if not continueRoutine: 
            # if we abort early the non-slip timer needs reset:
            routineTimer.reset() 
            break
        # will revert to True if at least one component still running
        continueRoutine = False  
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "status") and \
                    thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
            print np.shape(behResp)
            if et:
                elEndRec(el)
            core.quit()
        
        # refresh the screen
        # don't flip if this routine is over or we'll get a blank screen
        if continueRoutine:  
            win.flip()
        else: # this Routine was not non-slip safe so reset non-slip timer
            routineTimer.reset()
    
    #-------Ending Routine "trial"-------
    for thisComponent in trialComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)

if et:
    # File transfer and cleanup!
    pl.endRealTimeMode()
    el.setOfflineMode()						  
    pl.msecDelay(600) 

    #Close the file and transfer it to Display PC
    el.closeDataFile()
    el.receiveDataFile(edfFileName, edfFileName)
    os.rename(edfFileName, filePath + os.sep + edfFileName)
    el.close()

print "finished the experiment"

win.close()
core.quit()
