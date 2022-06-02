try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing
import math

# --------------------------- Define constants -------------------------

c = 299792458.0
lambda_0 = 1.0e-6
x_min = 0.0 * lambda_0
x_max = 120.0 * lambda_0
y_min = -20.0 * lambda_0
y_max = 20.0 * lambda_0
nx = 18000
ny = 1200
dx = (x_max - x_min) / nx
dy = (y_max - y_min) / ny
dt_multiplier = 0.99
dt = dt_multiplier * dx * dy / math.sqrt(dx**2 + dy**2) / c
T = lambda_0 / c
t_max = 120.0 * T
iterations = int(t_max / dt)
output_freq_0 = 1.0 * T / dt
output_freq_2 = 0.05 * T / dt
dump_freqs_0 = []
dump_freqs_2 = []
for i in range(iterations):
    dump_freqs_0.append(int(i * output_freq_0))
    dump_freqs_2.append(int(i * output_freq_2))
inputs = ['Grid']
update_freq = 1
root_dir = '/scratch/temp/valenpe7/single_cycle_2/2d/50'

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:

      data = coprocessor.CreateProducer(datadescription, 'Grid')

      #grid_writer = servermanager.writers.XMLPImageDataWriter(Input=grid, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
      #coprocessor.RegisterWriter(grid_writer, filename='./grid_%t.pvti', freq=1)

      #resample = ResampleToImage(Input=grid)
      #resample.SamplingDimensions = [int(nx / 3.0), int(ny / 6.0), int(nz / 6.0)]
      #resample.UseInputBounds = True
      #resample_writer = servermanager.writers.XMLPImageDataWriter(Input=resample, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
      #coprocessor.RegisterWriter(resample_writer, filename='resample_%t.pvti', freq=output_freq_1, paddingamount=5)

      #axis = PlotOverLine(Input=grid, Source='Line')
      #axis.Source.Point1 = [x_min, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0]
      #axis.Source.Point2 = [x_max, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0]
      #axis.Source.Resolution = nx
      #axis_writer = servermanager.writers.XMLPPolyDataWriter(Input=axis, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
      #coprocessor.RegisterWriter(axis_writer, filename='axis_%t.pvtp', freq=output_freq_2, paddingamount=5)

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

    def ProcessTriggers(self, datadescription):

      if datadescription.GetTimeStep() in dump_freqs_0:
        time = datadescription.GetTime()
        data = self.Pipeline.data
        data.UpdatePipeline(time)
        data.UpdatePipelineInformation()

        ex = PassArrays(Input=data)
        ex.PointDataArrays=['Electric Field - E_x (V/m)']
        ex.FieldDataArrays=['TimeValue']
        ex_writer = servermanager.writers.XMLPImageDataWriter(Input=ex, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
        ex_writer.FileName = 'ex_%d.pvti' % dump_freqs_0.index(datadescription.GetTimeStep())
        ex_writer.UpdatePipeline(time)

        ez = PassArrays(Input=data)
        ez.PointDataArrays=['Electric Field - E_z (V/m)']
        ez.FieldDataArrays=['TimeValue']
        ez_writer = servermanager.writers.XMLPImageDataWriter(Input=ez, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
        ez_writer.FileName = 'ez_%d.pvti' % dump_freqs_0.index(datadescription.GetTimeStep())
        ez_writer.UpdatePipeline(time)

        ne_back = PassArrays(Input=data)
        ne_back.PointDataArrays=['Number Density - electron (m^-3)']
        ne_back.FieldDataArrays=['TimeValue']
        ne_back_writer = servermanager.writers.XMLPImageDataWriter(Input=ne_back, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
        ne_back_writer.FileName = 'ne_back_%d.pvti' % dump_freqs_0.index(datadescription.GetTimeStep())
        ne_back_writer.UpdatePipeline(time)

        #ne_inj = PassArrays(Input=data)
        #ne_inj.PointDataArrays=['Number Density - electron_inj (m^-3)']
        #ne_inj.FieldDataArrays=['TimeValue']
        #ne_inj_writer = servermanager.writers.XMLPImageDataWriter(Input=ne_inj, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
        #ne_inj_writer.FileName = 'ne_inj_%d.pvti' % dump_freqs_0.index(datadescription.GetTimeStep())
        #ne_inj_writer.UpdatePipeline(time)

      if datadescription.GetTimeStep() in dump_freqs_2:
        time = datadescription.GetTime()
        data = self.Pipeline.data
        data.UpdatePipeline(time)
        data.UpdatePipelineInformation()
        '''
        axis = PlotOverLine(Input=data, Source='Line')
        axis.Source.Point1 = [x_min, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0]
        axis.Source.Point2 = [x_max, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0]
        axis.Source.Resolution = nx
        axis_writer = servermanager.writers.XMLPPolyDataWriter(Input=axis, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
        axis_writer.FileName = 'axis_%d.pvtp' % dump_freqs_2.index(datadescription.GetTimeStep())
        axis_writer.UpdatePipeline(time)
        '''
        axis = ExtractSubset(Input=data)
        axis.VOI = [1, 18001, 600, 600, 0, 0]
        axis.SampleRateI = 1
        axis.SampleRateJ = 1
        axis.SampleRateK = 1
        axis.IncludeBoundary = 0
        axis_writer = servermanager.writers.XMLPImageDataWriter(Input=axis, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
        axis_writer.FileName = 'axis_%d.pvti' % dump_freqs_2.index(datadescription.GetTimeStep())
        axis_writer.UpdatePipeline(time)

      '''
      if datadescription.GetTimeStep() in dump_freqs_1:
        time = datadescription.GetTime()
        data = self.Pipeline.data
        data.UpdatePipeline(time)
        data.UpdatePipelineInformation()
        resample = ResampleToImage(Input=data)
        resample.UseInputBounds = False
        resample.SamplingBounds = [x_min, x_max, -15.0e-6, +15.0e-6, 0.0, 0.0]
        resample.SamplingDimensions = [int(nx / 3.0), int(ny / 3.0), 1]
        writer = servermanager.writers.XMLPImageDataWriter(Input=resample, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
        writer.FileName = 'resample_%d.pvti' % dump_freqs_1.index(datadescription.GetTimeStep())
        writer.UpdatePipeline(time)

      if datadescription.GetTimeStep() in dump_freqs_2:
        time = datadescription.GetTime()
        data = self.Pipeline.data
        data.UpdatePipeline(time)
        data.UpdatePipelineInformation()
        axis = PlotOverLine(Input=data, Source='Line')
        axis.Source.Point1 = [x_min, (y_min + y_max) / 2.0, 0.0]
        axis.Source.Point2 = [x_max, (y_min + y_max) / 2.0, 0.0]
        axis.Source.Resolution = nx
        writer = servermanager.writers.XMLPPolyDataWriter(Input=axis, DataMode='Appended', EncodeAppendedData=0, HeaderType='UInt32', CompressorType='ZLib')
        writer.FileName = 'axis_%d.pvtp' % dump_freqs_2.index(datadescription.GetTimeStep())
        writer.UpdatePipeline(time)
      '''

  coprocessor = CoProcessor()
  freqs = {}
  for name in inputs:
    freqs[name] = [update_freq]
  coprocessor.SetUpdateFrequencies(freqs)
  coprocessor.SetRootDirectory(root_dir)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    coprocessor.ProcessTriggers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=True, image_quality=0, padding_amount=0)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "login3", 22222)
