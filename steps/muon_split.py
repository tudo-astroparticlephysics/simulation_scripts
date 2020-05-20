### Define functions for muon splitting ###


def muon_splitter(tray, name):
    from icecube import dataclasses, icetray

    def dir_pos(frame):
        frame.Put('NewDirection', dataclasses.I3Direction())
        frame.Put('NewPosition', dataclasses.I3Position())		
    
    def new_tree(frame):
        frame.Put('NewTree', dataclasses.I3MCTree())


    tray.AddModule(dir_pos, 'addDirPos', Streams=[icetray.I3Frame.DAQ])
    tray.AddModule(new_tree, 'addTree', Streams=[icetray.I3Frame.DAQ])
	

   
