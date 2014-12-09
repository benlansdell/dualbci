#Import network
network import file indexColumnTargetInteraction=1 indexColumnSourceInteraction=2 file="/home/lansdell/projects/bci/matlab/worksheets/12_3_2014/GLMGrangerB.xml.graphml"

#Import and set style
vizmap load file file="/home/lansdell/projects/bci/matlab/worksheets/12_3_2014/GrangerStyle.xml"
vizmap apply styles=DirectedGranger

#Set layout
layout attributes-layout NodeAttribute=cluster maxwidth=400
#Set view to fit display
view fit content
#Save 
view export OutputFile="/home/lansdell/projects/bci/matlab/worksheets/12_3_2014/ClusteredGrangerB.pdf" options=PDF

#Set layout
layout attribute-circle
#Set view to fit display
view fit content
#Save 
view export OutputFile="/home/lansdell/projects/bci/matlab/worksheets/12_3_2014/CircleGrangerB.pdf" options=PDF

#Quit
#command quit
