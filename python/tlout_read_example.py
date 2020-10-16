import tlout_read as tlor

def callback(timelines, shift, maxed_out_time, extinction_time):
    print(shift, maxed_out_time, extinction_time,timelines)

tlor.read("../l2_mu3_pinf8.tlout", callback)
