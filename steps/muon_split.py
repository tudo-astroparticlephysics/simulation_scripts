### Define functions for muon splitting ###


def muon_splitter(tray, name):
    from icecube import dataclasses, icetray

    def dir_pos_tree(frame):
        frame.Put('NewDirection', dataclasses.I3Direction())
        frame.Put('NewPosition', dataclasses.I3Position())		
        frame.Put('NewTree', dataclasses.I3MCTree())


    tray.AddModule(dir_pos_tree, 'addDirPos', Streams=[icetray.I3Frame.DAQ])
	

   
