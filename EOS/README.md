# How to use?

```
import OBS_rotate
EOSROT = OBS_rotate(OriDat="path_to_obs_ori_dat)")
rotated_stream = EOSROT.rotate_EOS_station(st=unrotated_stream)
```

unrotated_stream : original seismogram trace from CWB with channel code=1,2,3

rotated_stream : rotated seismogram trace with channel code=Z,N,E