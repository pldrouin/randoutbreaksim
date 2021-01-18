import tlout_read as tlor

def callback(timelines, shift, maxed_out_time, extinction_time):
    print(shift, maxed_out_time, extinction_time,timelines)

tlor.read("no_ct_no_mask_1M_paths_30day_tmax.tlout", callback)
