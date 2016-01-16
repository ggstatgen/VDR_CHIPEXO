
def plot_density(data_points, outfile, title = ""):
    pass

def plot_datalines(data_points, outfile, xlabels = None, title = "", xlab = '', ylab = ''):
    graphics = importr('graphics')
    grdevices = importr('grDevices')
    r_base = importr('base')

    grdevices.pdf(outfile)
    

    graphics.plot(data_points[0], data_points[1], type="o", xlab = xlab, ylab= ylab, main=title)
    
    if xlabels != None:
        xlab = StrVector(xlabels)

        graphics.axis(1, at = IntVector(range(1,len(xlabels)+1)), lab=xlab)
        graphics.axis(2)

    grdevices.dev_off()

    return



def plot_barplot(data_points, outfile, xlabel = '', ylabel = '', title = "", line_width=2):
    graphics = importr('graphics')
    grdevices = importr('grDevices')
    r_base = importr('base')


    grdevices.postscript(outfile)

    X = FloatVector(data_points[0])
    Y = FloatVector(data_points[1])

    graphics.plot(X, Y, xlab=xlabel, ylab=ylabel, type='n', main = title)
    graphics.segments(X, Y, X,FloatVector([0]*len(data_points[0])), lwd=line_width)
    grdevices.dev_off()

    return


def plot_multiple_datalines(data_points,legend_col, outfile, xlabels = None, title = "", xlabel = '', ylabel =''):
    graphics = importr('graphics')
    grdevices = importr('grDevices')
    r_base = importr('base')


    grdevices.pdf(outfile)

   # graphics.par(new=True)
    for i in range(len(data_points)):
        
        data = data_points[i]
        graphics.plot(data[0], data[1], axes = False, type="o", col = legend_col[1][i], xlab = xlabel, ylab = ylabel, main = title)
        graphics.par(new=True)


    
    print legend_col[1]
    graphics.legend("topright", leg=StrVector(legend_col[0]), fill=StrVector(legend_col[1]), cex=0.5, bty="n")

    if xlabels != None:
        xlab = StrVector(xlabels)

        graphics.axis(1, at = IntVector(range(1,len(xlabels)+1)), lab=xlab)
        graphics.axis(2)

    grdevices.dev_off()


def plot_histogram(data_points, outfile, breaks, title = "", xlabel = '', ps = False, stats_bool = True, log_scale = ""):

    import tempfile 
    from subproc import run
                                               
    import os
    #cwd = os.getcwd()
    scripts = '''
X <- read.table("tmp.txt")
postscript("%s")
hist(X[,],breaks=%s, main="%s", xlab="%s", col='red')
dev.off()
'''  % (outfile, str(breaks), title, xlabel)

    #print scripts
    mTempdir = tempfile.mkdtemp()
    #print data_points
    input_file = file("%s/tmp.txt"%mTempdir,'w')
    input_file.write("\n".join([str(x) for x in data_points]))
    input_file.close()
    script_file = file("%s/script.R"%mTempdir,'w')
    script_file.write(scripts)
    script_file.close()


    print run("cd %s;Rscript script.R"%mTempdir)
   # print mTempdir
    run("rm -rf %s" % mTempdir) 


def plot_histogram_fancy(scripts, data_points):
    import tempfile
    from subproc import run

    import os
    #cwd = os.getcwd()                                                                                                                            

    #print scripts                                                                                                                                 
    mTempdir = tempfile.mkdtemp()
    #print data_points                                                                                                                            
 
    input_file = file("%s/tmp.txt"%mTempdir,'w')
    input_file.write("\n".join([str(x) for x in data_points]))
    input_file.close()
    script_file = file("%s/script.R"%mTempdir,'w')
    script_file.write(scripts)
    script_file.close()

    print run("cd %s;Rscript script.R"%mTempdir)# print mTempdir                                                                                
    run("rm -rf %s" % mTempdir)



if __name__ == "__main__":
    import random
    plot_histogram([random.randint(0,5) for i in range(1000)], "testplot.pdf", (0,5,1),'density of random integers')
