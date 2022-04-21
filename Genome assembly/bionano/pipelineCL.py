import sys
import atexit
import os
import utilities as util
import traceback

"""
@package pipelineCL 
Entry point of Pipeline for command line

Usage: python pipelineCL.py -h
"""

    
#this is moved here from utilities/Pipeline.py because those files are used outside of the pipeline
@atexit.register
def on_exit():
	try:
            util.LogStatus("pipeline", "exit", "%d" % os.getpid())
	except:
            pass


if __name__ == "__main__":
    import Pipeline
    varsP = Pipeline.varsPipeline()

#    print "Starting _main_ try:" # DEBUG DEBUG
#    sys.stdout.flush()

    try :
            #    print "Before prerunChecks: varsP.autoNoise=", varsP.autoNoise # DEBUG
            varsP.prerunChecks()
    
            #    print "After prerunChecks: varsP.autoNoise=", varsP.autoNoise # DEBUG
            #    varsP.updatePipeReport("After prerunChecks: varsP.autoNoise= %d\n\n" % varsP.autoNoise) # DEBUG

            print('  Prerun Tests:\n\t%d ERRORS\n\t%d WARNINGS\n' % (varsP.error, varsP.warning))
            if varsP.error or varsP.warning:
                    #print(varsP.message)
                    varsP.printMessage()
            if varsP.error:
                    print('  EXITING: See errors') 
                    sys.exit(1)
    
            #    print "Before calling DNPipeline: varsP.autoNoise=", varsP.autoNoise # DEBUG
            #    varsP.updatePipeReport("Before calling DNPipeline: varsP.autoNoise= %d\n\n" % varsP.autoNoise) # DEBUG

            dnpipeline = Pipeline.DNPipeline()

            #    print "Before calling DNPipeline.run(varsP): varsP.autoNoise=", varsP.autoNoise # DEBUG
            #    varsP.updatePipeReport("Before calling DNPipeline.run(varsP): varsP.autoNoise= %d\n\n" % varsP.autoNoise) # DEBUG
            os.putenv('OMP_WAIT_POLICY', 'PASSIVE')
            dnpipeline.run(varsP)

    except Exception, e:
            print "Python Exception : %s" % str(e)
            print Exception

            print "exception stack trace:"
            traceback.print_exc()

            sys.stdout.flush()
            sys.stderr.flush()

            os._exit(1) # terminates all threads
            sys.exit(1)
