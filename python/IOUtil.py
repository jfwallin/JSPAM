
class IOUtil:

  @staticmethod
  def outputParticles(filename,x0):
    fo = open(filename,'w')
    IOUtil.outputParticlesToFile(fo,x0)


  @staticmethod
  def outputParticlesToFile(fo,x0):
    size = len(x0)
    for i in range(size):
      dtmp = x0[i]
      for j in range(6):
        fo.write(IOUtil.formatDouble(dtmp[j]))
      # new line 
      fo.write("\n")
  
    fo.close()

  @staticmethod
  def formatDouble(num):
    return "%16.8e"%(num)
