from obspy.core import Stream
import numpy as np
import numpy.linalg as LA


class OBS_rotate():
    """
    This class rotate CWB EOS stations from input stream and return a new stream stored rotated station with new channel id.
    
    Program will automatic trim and fill_value for input stream to rotate.

    Make sure input EOS traces with the same starttime and endtime.
    """
    def __init__(self, OriDat="obs_ori.dat"):
        # path to orientation info file
        self.OriDat = OriDat
    
    def data_pre_process(self, st):
        tmp_st = st.select(network='TW', location="20").copy()
        # check EOS station
        if len(tmp_st) == 0:
            print(f"Can't find any location code=20 in TW network, no EOS stations!")
            return None
        tmp_st.merge(fill_value=0)
        sttime = []
        edtime = []
        for tr in st:
            sttime.append(tr.stats.starttime)
            edtime.append(tr.stats.endtime)
        sttime.sort()
        edtime.sort()
        # trim lastest starttime and earliest endtime
        trim_sttime = sttime[-1]
        trim_edtime = edtime[0]
        tmp_st.trim(trim_sttime, trim_edtime)
        return tmp_st
    
    def rotate_EOS_station(self, st):
        # set all possible channel pair in EOS2~EOSA
        CHAN_PAIR = (['HL1','HL2','HL3'], ['HN1','HN2','HN3'], ['EH1','EH2','EH3'], ['HA1','HA2','HA3'])

        # pre_process
        tmp_st = self.data_pre_process(st)
        if tmp_st == None:
            return None

        # find EOS station in time period
        self.st_time = tmp_st[0].stats.starttime
        self.get_station_info()

        # find EOS stations and rotate it
        new_EOS_st = Stream()
        for STN in self.STNs:
            for PAIR in CHAN_PAIR:
                EOS_st = Stream()
                for CHN in PAIR:
                    EOS_st += tmp_st.select(station=STN).select(channel=CHN)
                
                if len(EOS_st) == 3:
                    PHI = np.deg2rad(self.STNs_Roll[STN])
                    PSI = np.deg2rad(self.STNs_Pitch[STN])
                    azimuth = self.STNs_Az[STN] - 90
                    self.rotate_seis(EOS_st, PHI, PSI, azimuth)
                    # change channel code from 123 to ZNE
                    tr_id = EOS_st[0].stats.channel[0:2]
                    EOS_st[0].stats.channel = f"{tr_id}Z"
                    EOS_st[1].stats.channel = f"{tr_id}N"
                    EOS_st[2].stats.channel = f"{tr_id}E"
                    # combine to output stream
                    new_EOS_st += EOS_st
        return new_EOS_st
    
    def make_YAW_PITCH_ROLL(self, x_vector, y_vector, z_vector):
        """
        calculate rotate degree from input x, y, z vector (mean of daily seismogram).
        """
        YZ_len  = LA.norm(np.array([y_vector, z_vector]))
        PHI = np.arctan2(y_vector, z_vector)
        PSI = np.arctan2(-x_vector, YZ_len)
        return PHI, PSI
    
    def make_rotate_matrix(self, PHI, PSI, azimuth):
        """
        pre-multiply three rotation matrix for rotate seismogram
        """
        COS_PHI = np.cos(-PHI)
        SIN_PHI = np.sin(-PHI)
        
        ROTATEM_PHI = np.array([[1.0, 0.0, 0.0], \
                            [0.0, COS_PHI,  SIN_PHI], \
                            [0.0, -SIN_PHI,  COS_PHI]])
        
        COS_PSI = np.cos(-PSI)
        SIN_PSI = np.sin(-PSI)
        
        ROTATEM_PSI = np.array([[ COS_PSI, 0.0, -SIN_PSI], \
                                [ 0.0, 1.0, 0.0], \
                                [ SIN_PSI, 0.0, COS_PSI]])
        
        
        COS_AZI = np.cos(np.deg2rad(azimuth))
        SIN_AZI = np.sin(np.deg2rad(azimuth))
        
        ROTATEM_azimuth = np.array([[ COS_AZI, SIN_AZI, 0.0], \
                                    [-SIN_AZI, COS_AZI, 0.0], \
                                    [     0.0,     0.0, 1.0]])
        
        ROTATE_M_TMP = np.matmul(ROTATEM_PSI, ROTATEM_PHI)
        ROTATE_M = np.matmul(ROTATEM_azimuth, ROTATE_M_TMP)
        return ROTATE_M
    
    def rotate_seis(self, st, PHI, PSI, azimuth):
        """
        rotate seismogram
        """
        tmp = np.append([st[2].data.data, st[1].data.data], [st[0].data.data], axis=0)
        SEIS_3D = tmp.reshape(3, st[2].stats.npts)
        ROTM_3D = self.make_rotate_matrix(PHI, PSI, azimuth)

        SEIS_NEW = np.matmul(ROTM_3D, SEIS_3D)

        st[2].data = SEIS_NEW[0]
        st[1].data = SEIS_NEW[1]
        st[0].data = SEIS_NEW[2]
        
        return
        
    def get_station_info(self):
        """
        get station orientation info from input file.
        """
        ORI_DAT = open(self.OriDat, 'r')
        TIME = self.st_time.year*10000+self.st_time.month*100+self.st_time.day
        self.STNs = []
        self.STNs_Roll   = {}
        self.STNs_Pitch  = {}
        self.STNs_Az     = {}
        
        for line in ORI_DAT:
            STN = line[0:5].strip()
            period_start = int(line[37:45])
            period_stop  = int(line[46:54])
            
            if STN not in self.STNs:
                if TIME >= period_start and TIME <= period_stop:
                    self.STNs.append(STN)
                    self.STNs_Roll[STN]  = float(line[8:18])
                    self.STNs_Pitch[STN] = float(line[19:29])
                    self.STNs_Az[STN]    = float(line[30:36])
        ORI_DAT.close()
        return
