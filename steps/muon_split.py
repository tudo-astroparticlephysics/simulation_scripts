### Define functions for muon splitting ###



def angle_pos(frame, angle, pos):
    from icecube import dataclasses

    frame.Put('ChangeAngleZe', dataclasses.I3Double(angle[0]))
    frame.Put('ChangeAngleAz', dataclasses.I3Double(angle[1]))
    frame.Put('ChangePosX', dataclasses.I3Double(pos[0]))
    frame.Put('ChangePosY', dataclasses.I3Double(pos[1])) 	
    frame.Put('ChangePosZ', dataclasses.I3Double(pos[2])) 
    print(frame)
