'''
This is the python version of pit removal, credit to Stephen Jackson
'''

import numpy as np
import os
import gdal
import matplotlib.pyplot as plt
from qgis.core import QgsApplication, QgsMapLayer, QgsProject, QgsRasterLayer



class _PitRemoval(object):
    '''
    Inputs
    -------------
    dem_map: str, file path,  DEM raster data
    layer: QGIS raster layer object
    mode: str, "CUT", "BAL", "MIN_COST", "MIN_CELL",
         ('cut' for cut only, 'bal' for minimize sum(delta z),
         'mincost' for minimize sum(abs(delta z)),
         'mincell' for minimize number of cells modified)
    step_size: float, vertical step size
    dem_dst: str, destination for produced DEM file

    Attributes:
    -------------
    geo_info: numpy structured data
            dem: dem value for given basin
            border: border map
            flow_dir: flow direction defined from East at 0 and northwest at 128
            flooded: a cell is flooded or not with 0 unflooded, 1 flooded 2 has confirmed descending path to the outlet
            checked: boolean value for temporary storage.

    main_queue: Queue Object
            id: number of cell to be processed
            elevation: elevation extracted from dem

    depression: List, depression extent from a given ID

    dry_neighbors:Queue Object

    Methods:
    --------------
    _init_main_queue: firstly run this to initialize main queue
    .
    itermate_main_queue: main function to execute OptimizedPitRemoval

    Return:
    --------------
    None
    '''
    def __init__(self,layer= None,  dem_map=None, mode="MIN_COST", step_size=0.1, dem_dst=None):
        self.dem_map= dem_map
        if layer is not None:
            if not isinstance(layer, QgsRasterLayer):
                raise ValueError("Invalid input type, expected QgsRasterLayer, but get %s"%str(type(layer)))
            else:
                self.dem= self._read_dem_layer(layer)

        if dem_map is not None:
            if not (isinstance(dem_map, str) or os.path.exists(dem_map)):
                raise ValueError('File doesnot exist!')
            else:
                self.dem= self._read_dem(self.dem_map) # dem data
            # print(' load in raster map: \ncrs:  %s\nheight:  %d\nwidth:  %d\ncount:  %d\ntransform:\n%s'%(self.crs,
            #  self.width, self.height, self.count, self.transform))
        if (dem_map is None and layer is None):
            raise ValueError('at least provide one data source!')

        self.mode= mode
        self.step_size= step_size
        self.dem_dst= dem_dst
        self.num_rows, self.num_cols= self.dem.shape
        self.tot_size= self.num_rows*self.num_cols
        self._init_queue(self.dem)

    def iterate_main_queue(self):
        '''Main function to execute OptimizedPitRemoval'''
        num_pit=0
        ConfirmDescend= False

        while(not self.main_queue.empty()): # executes until emptying all queues
            cur_cell= self.main_queue.top()
            self.main_queue.pop()
            # check if it is local minimum
            if self.isLocalMinimum(cur_cell.id):
                if self.mode=='CUT':
                    self.cut_to_elevation(cur_queue.id)
                else:
                    self.pit_removal_hybrid(cur_queue.id)

                num_pit+=1
            else:
                m_cur, n_cur= self._id_transform(cur_queue.id)
                if self.geo_info['flooded'][m_cur, n_cur]==1:
                    neighbors= self._get_neighbors(cur_cell.id)
                    for i in range(8):
                        if neighbors[i]>-1:
                            m_nei, n_nei= self._id_transform(neighbors[i])
                            if self.geo_info['flooded'][m_nei, n_nei]==2 and self.geo_info['dem'][m_nei, n_nei]<= cur_cell.elev:
                                self.geo_info['flooded'][m_nei, n_nei]=2
                                break

            #Add unflooded neighbors to Main Queue and identify the direction the flooding came from
            self.get_dry_neighbors(cur_cell.id)
            while not self.neighbor_queue.empty():
                cur_neighbor= self.neighbor_queue.top()
                self.neighbor_queue.pop()

                ConfirmDescend= False
                m_nei, n_nei= self._id_transform(cur_neighbor.id)
                m_cur, n_cur= self._id_transform(cur_cell.id)
                if cur_neighbor.elev>= cur_cell.elev and self.geo_info['flooded'][m_cur, n_cur]==2:
                    ConfirmDescend= True
                self._add_main_queue(cur_neighbor.id, ConfirmDescend)
                direction= self._flow_dir(cur_neighbor.id, cur_cell.id)
                self.geo_info['flow_dir'][m_nei, n_nei]= direction

    def pit_removal_hybrid(self, ID):
        '''
        Description:
        ----------------------
        This handles the case where mode= BAL or MIN_COST
        Cost is definede as the difference between the original terrain for each point summed accross all modified points.
        BAL minimizes the sum of cost across all cells modified for each depression (i.e. tries to get Cut Volume=Fill Volume)
        MIN_COST minimizes the sum of the absolute value of cost across all cells modified for each depression (i.e. tries to have the least disturbance to the original terrain)
        BAL will *always* have a mix of cut and fill (for suitably small step size)
        MIN_COST will *sometimes* have a mix of cut and fill, but some pits may be best removed using only cut or only fill
        Filling modifies a 2-D region (depression). Cutting modifies a 1-D path to an outlet. Thus MIN_COST will often result in more cut than fill.


        '''

        m,n= ID//self.num_cols, ID%self.num_rows
        pit_elev= self.geo_info['dem'][m,n]

        crest_elev= self.get_crest_elev(ID)
        # get depression extents and set pit elevation
        self.get_depression_extent(ID, crest_elev)
        # create Cut function
        cut= self.create_cut_func(ID, crest_elev)
        # create fill function
        fill= self.create_fill_func(ID, crest_elev)

        if pit_elev<crest_elev:
            # find desired fill elevation
            cut, fill, cost, ideal_fill_level= self.get_ideal_fill_level(ID, crest_elev, cut, fill)

            self.fill_to_elevation(ID, ideal_fill_level)

            self.cut_to_elevation(ID)



    def get_crest_elev(self, ID):
        m,n= self._id_transform(ID)
        reached_outlet=0
        crest_elev= self.geo_info['dem'][m,n]

        while not reached_outlet:
            nextID= self.trace_flow(ID, self.geo_info['flow_dir'][m,n])
            next_m, next_n= self._id_transform(nextID)
            if NextID<0:
                reached_outlet=1
            elif np.isnan(self.geo_info['dem'][next_m, next_n]):
                reached_outlet=1
            elif self.geo_info['dem'][next_m, next_n]< crest_elev and self.geo_info['flooded'][next_n, next_m]==2:#NextID is lower than Pit and NextID on confirmed descending path to outlet
                reached_outlet=1
            else:
                if self.geo_info['dem'][next_m, next_n]> crest_elev and not np.isnan(self.geo_info['dem'][next_m, next_n]):
                    crest_elev= self.geo_info['dem'][next_m, next_n]
            ID= next_ID

        return elevation

    def get_depression_extent(self, ID: int, elev_crest: float):
        # To make an elevation ordered list of every cell in the depression
        #A compound depression (neighboring pit with separating ridge lower than crest elevation) is treated as a separate depression. That pit will be removed later.
        depression_queue= Queue(name='depression')
        self.depression= []
        cur_pnt_m, cur_pnt_n= self._id_transform(ID)
        cur_pnt= Point(ID, self.geo_info['dem'][cur_pnt_m, cur_pnt_n])
        depression_queue.push(cur_pnt)
        self.depression.append(cur_pnt_id)

        while (not depression_queue.empty()):
            cur_pnt= depression_queue.top()
            depression_queue.pop()
            cur_id= cur_pnt.id

            for i in range(8):
                check, neigh_pnt= self.check_cell(cur_id, i, 0, crest_elev)
                if check:
                    depression_queue.push(neigh_pnt)
                    self.depression.append(neigh_pnt.id)

        self.geo_info['checked']= np.zeros((self.num_rows, self.num_cols), dtype=bool)

    def create_cut_func(self, ID: int, crest_elev: float):
        '''
        Desciption:
        -------------------
        Create a list of elevations from the Pit to the Crest, incremented by step size
        For each elevation in list, sum all positive difference between Terrain elev and List Elev for all cells along path from pit to outlet
        CutFunction represents the (Cut Volume / Cell Area)

        Return:
        -------------------
        cut: dict, cut values for each cell along the descending path of given ID
        '''
        self.cut= dict()
        reached_outlet= False

        m,n= self._id_transform(ID)
        cell_elev= self.geo_info['dem'][m,n]
        pit_elev= self.geo_info['dem'][m,n]
        cur_id= ID

        cur_step= self.geo_info['dem'][m,n]
        while cur_step< crest_elev:
            self.cut[cur_step]= 0
            cur_step+= self.step_size

        self.cut[crest_elev]= 0

        while (not reached_outlet):
            m_cur,n_cur= self._id_transform(cur_id)
            cur_flow_dir= self.geo_info['flow_dir'][m,n]
            next_id= self.trace_flow(cur_id, cur_flow_dir)
            m_next, n_next= self._id_transform(next_id)
            if next_id<0: # if is border cell
                reached_outlet= 1
            elif np.isnan(self.geo_info['dem'][m_cur,n_cur]):
                reached_outlet= 1
            elif (self.geo_info['dem'][m_cur,n_cur]< pit_elev) and self.geo_info['flooded'][m_next, n_next]:
                reached_outlet= 1
            else:
                cell_elev= self.geo_info['dem'][m_next, n_next]
                cur_step= pit_elev

                while cur_step< cell_elev:
                    old_cut= [cut[cur_step] if len(cut[cur_step]) else 0][0]
                    cut[cur_step]= old_cut+ cell_elev- cur_step
                    cur_step+= self.step_size

            cur_id= next_id

        return cut

        def create_fill_func(self, ID: int, crest_elev: float):
            '''
            Description:
            ------------------
            Create a list of elevations from the Pit to the Crest, incremented by step size
            For each elevation in list, sum all positive difference between List elev and Terrain Elev for all cells within the depression
            FillFunction represents the (Fill Volume / Cell Area)

            Return:
            ------------------
            fill: dict, fill cost of descending path for given ID
            '''
            m_cur, n_cur= self._id_transform(ID)
            pit_elev= self.geo_info['dem'][m_cur, n_cur]
            cur_step= pit_elev
            fill= dict()
            while cur_step< crest_elev:
                fill[cir_step]= 0
                cur_step+= self.step_size

            fill[crest_elev]= 0
            # calculate cost
            for i in range(len(self.depression)):
                cur_id= self.depression[i]
                m_cur, n_cur= self._id_transform[cur_id]
                cur_ground_elev= self.geo_info['dem'][m_cur, n_cur]

                cur_step= pit_elev
                while cur_step< crest_elev:
                    if cur_step> cur_ground_elev:
                        old_fill= fill[cur_step]
                        fill[cur_step]= old_fill+ cur_step -cur_ground_elev

                    cur_step+= self.step_size

                    if cur_step> crest_elev:
                        old_fill= fill[crest_elev]
                        fill[crest_elev]= old_fill+ crest_elev- cur_ground_elev

            return fill

    def get_ideal_fill_level(self,ID: int, crest_elev: float, cut: dict, fill: dict):
        '''
        Description:
        -----------------
        MODE: 'BAL': balancing desired to find the minimum difference of fill and cut
              'MIN_COST': find the point of minimum cost

        Return:
        -----------------

        '''
        m,n= self._id_transform(ID)
        pit_elev= self.geo_info['dem'][m,n]
        def BAL(pit_elev, cut, fill):
            cur_step= pit_elev
            cur_cut_cost= cut[cur_step]
            cur_fill_cost= fill[cur_step]
            min_diff= abs(cur_fill_cost- cur_cut_cost)
            fill_level= cur_step

            #Optimize:
            while cur_step< crest_elev:
                cur_cut_cost= cut[cur_step]
                cur_fill_cost= fill[cur_step]
                cur_diff= abs(cur_fill_cost- cur_cut_cost)
                if cur_diff< min_diff:
                    min_diff= cur_diff
                    fill_level= cur_step

                cur_step+= step_size

            return cut, fill, min_diff, fill_level


        def MIN_COST(pit_elev, cut, fill):
            cur_step= pit_elev
            min_cost= cut[cur_step]
            fill_level= cur_step
            if min_cost>0:
                while cur_step<crest_elev:
                    cur_cut_cost= cut[cur_step]
                    cur_fill_cost= fill[cur_step]
                    cur_total_cost= cur_cut_cost+ cur_fill_cost

                    if cur_total_cost< min_cost:
                        min_cost= cur_total_cost
                        fill_level= cur_step

                    cur_step+= self.step_size

            return cut, fill, min_cost,fill_level

        switcher= {'BAL': BAL,
                    'MIN_COST': MIN_COST}

        cut, fill, cost, fill_level= switcher[self.mode](pit_elev, cut, fill)

        return cut, fill, cost, fill_level

    def fill_to_elevation(self, ID: int, fill_level: int):
        m,n= self._id_transform(ID)
        if self.geo_info['dem'][m,n]<fill_level and not np.isnan(self.geo_info['dem'][m,n]): self.geo_info['dem'][m,n]=fill_level

        for i in range(len(self.depression)):
            cur_id= self.depression[i]
            m_cur, n_cur= self._id_transform(cur_id)
            if self.geo_info['dem'][m_cur, n_cur]<fill_level and not np.isnan(self.geo_info['dem'][m_cur, n_cur]): self.geo_info['dem'][m_cur, n_cur]=fill_level

    def cut_to_elevation(self, ID: int):
        m,n= self._id_transform(ID)
        pit_elev= self.geo_info['dem'][m,n]
        reached_outlet=False

        cur_id= ID

        while (not reached_outlet):
            m_cur, n_cur= self._id_transform(cur_id)
            next_id= self.trace_flow(cur_id, self.geo_info['flow_dir'][m_cur, n_cur])
            m_next, n_next= self._id_transform(next_id)
            if next_id<0:
                reached_outlet=1
            elif np.isnan(self.geo_info['dem'][m_next, n_next]):
                reached_outlet=1
            elif self.geo_info['dem'][m_next, n_next]<pit_elev and self.geo_info['flooded'][m_next,n_next]==2:
                reached_outlet=1
            else:
                if self.geo_info['dem'][m_next, n_next]>pit_elev and not np.isnan(self.geo_info['dem'][m_next, n_next]):
                    self.geo_info['dem'][m_next, n_next]= pit_elev
                self.geo_info['flooded'][m_next, n_next]=2

            cur_id= next_id


    def check_cell(self, ID: int, direction: int, cur_neigh_id: int, crest_elev: float):

        def northwest(ID, num_cols):
            if not (ID<num_cols or ID %num_cols==0):
                cur_neigh_id= ID-1-num_cols
                check= True

            return cur_neigh_id, check

        def north(ID, num_cols):
            if not ID<num_cols:
                cur_neigh_id= ID- num_cols
                check= True

            return cur_neigh_id, check

        def northeast(ID, num_cols):
            if not (ID< num_cols or (ID+1)%num_cols==0):
                cur_neigh_id= ID+1-num_cols
                check= True

            return cur_neigh_id, check

        def east(ID, num_cols):
            if (ID+1) % num_cols:
                cur_neigh_id= ID+1
                check= True

            return cur_neigh_id, check

        def south(ID, num_cols):
            if not (self.tot_size- ID < num_cols+ 1):
                cur_neigh_id= ID + num_cols
                check = True

            return cur_neigh_id, check

        def southwest(ID, num_cols):
            if (not self.tot_size-ID< num_cols+1 and ID % num_cols):
                cur_neigh_id= ID-1+num_cols
                check= True

            return cur_neigh_id, check

        def west(ID, num_cols):
            if ID%num_cols:
                cur_neigh_id= ID-1
                check=True
            return cur_neigh_id, check

        check= False
        switcher= {0: northwest(ID, self.num_cols),
                   1: north(ID, self.num_cols),
                   2: northeast(ID, self.num_cols),
                   3: east(ID, self.num_cols),
                   4: southeast(ID, self.num_cols),
                   5: south(ID, self.num_cols),
                   6: southwest(ID, self.num_cols),
                   7: southwest(ID, self.num_cols)}
        m,n= self._id_transform(ID)
        cur_neigh_id, check= switcher[direction]
        cur_neigh_m, cur_neigh_n= self._id_transform(cur_neigh_id)
        cur_neigh_pnt= Point(ID, self.geo_info['dem'][cur_neigh_m, cur_neigh_n])
        if check:
            if checked[cur_neigh_id]!=0 and cur_neigh_pnt.elev< crest_elev and cur_neigh_pnt.elev> self.geo_inof['dem'][m,n]:
                check= False
            self.geo_info['checked'][cur_neigh_m, cur_neigh_n]=1

        return check, cur_neigh_pnt



    def trace_flow(self, ID_from: int, flow_direction: int) -> int:
        '''Inversed operation of calculating flow direction'''
        if(flow_direction == 1):
            ID_to = ID_from + 1;
        elif(flow_direction == 2):
            ID_to = ID_from + 1 + numCols;
        elif(flow_direction == 4):
            ID_to = ID_from + numCols;
        elif(flow_direction == 8):
            ID_to = ID_from - 1 + numCols;
        elif(flow_direction == 16):
            ID_to = ID_from - 1;
        elif(flow_direction == 32):
            ID_to = ID_from - 1 - numCols;
        elif(flow_direction == 64):
            ID_to = ID_from - numCols;
        elif(flow_direction == 128):
            ID_to = ID_from + 1 - numCols;
        else:
            ID_to = -1;

        return ID_to

    def get_dry_neighbors(self, ID: int):
        dry_neighbor= Point()
        neighbors= self._get_neighbors(ID)
        self.neighbor_queue= Queue()

        for i in range(8):
            if neighbors[i]>-1:
                m_nei, n_nei= self._id_transform(neighbors[i])
                if not self.geo_info['flooded'][m_nei, n_nei]:
                    dry_neighbor.id= Neighbors[i]
                    dry_neighbor.elev= self.geo_info['dem'][m_nei, n_nei]
                    self.neighbor_queue.push(dry_neighbor)

    def __repr__(self):

        return (' load in raster map: \ncrs:  %s\nheight:  %d\nwidth:  %d\ncount:  %d\ntransform:\n%s'%(self.crs,
             self.width, self.height, self.count, self.transform))


    def _init_queue(self, terrain):
        add_queue=0
        # initialize numpy object storing dem, border, flow_direction
        self.geo_info= {'dem': self.dem,
                        'border': np.zeros((self.num_rows, self.num_cols), dtype=np.int16),
                        'flow_dir': np.zeros((self.num_rows, self.num_cols), dtype=np.int16),
                        'flooded':np.zeros((self.num_rows, self.num_cols),dtype=np.int16),
                        'checked': np.zeros((self.num_rows, self.num_cols), dtype=bool)}
        self.main_queue= Queue(name='MainQueue')
        self.geo_info['dem']= self.dem
        # assign border to be 1, else 0
        self.geo_info['border'][0,:]=1; self.geo_info['border'][-1,:]=1; self.geo_info['border'][:,0]=1; self.geo_info['border'][:,-1]=1
        del self.dem
        # self.geo_info['dem'][2,2]=np.nan
        for i in range(self.tot_size):
            m= i//self.num_cols
            n= i% self.num_rows
            novalue= self._neighbor_wo_value(i)
            if self.geo_info['border'][m,n]:
                # print('border')
                add_queue= 1
            elif novalue:
                # print('no value')
                add_queue=1
            if add_queue:
                self._add_main_queue(i, True) #update main_queue
                self.geo_info['flow_dir'][m,n]= 0
            add_queue= 0

        return None


    def _add_main_queue(self,ID: int,ConfirmDescend: bool):
        m= ID//self.num_cols
        n= ID % self.num_rows
        if self.geo_info['flooded'][m,n]==0:
            self.main_queue.id.append(ID)
            self.main_queue.elev.append(self.geo_info['dem'][m,n])
            if ConfirmDescend:
                self.geo_info['flooded'][m,n]=2
            else:
                self.geo_info['flooded'][m,n]=1



    def _neighbor_wo_value(self,ID: int) -> int:
        novalue=0
        # check whether its neighbor has no value
        neighbors= self._get_neighbors(ID)
        for i in range(len(neighbors)):
            m,n= int(neighbors[i]//self.num_cols), int(neighbors[i] %self.num_rows)
            # print(m,n)
            if neighbors[i]==-1:
                novalue=1
            elif np.isnan(self.geo_info['dem'][m,n]):
                novalue=1
                self.geo_info['flow_dir'][m,n]= self._flow_dir(ID, neighbors[i])

        return novalue

    def _get_neighbors(self, ID):

        # Neighbors is a 0-7 vector, defined clockwise from Northwest
        # Returns the ID value for the eight neighbors, with -1 for cell off the grid
        Neighbors= np.zeros(8)
        #Northwest
        if ID<self.num_cols:
            Neighbors[0]= -1
        elif (ID % self.num_cols==0):
            Neighbors[0]=-1
        else:
            Neighbors[0]= ID-1-self.num_cols

        #North
        if ID<self.num_cols:
            Neighbors[1]=-1
        else:
            Neighbors[1]= ID-self.num_cols

        # Northeast
        if (ID<self.num_cols):
            Neighbors[2]= -1
        elif ((ID+1) % self.num_cols==0):
            Neighbors[2]= -1
        else:
            Neighbors[2]= ID+1-self.num_cols

        # East
        if ((ID+1)%self.num_cols==0):
            Neighbors[3]= -1
        else:
            Neighbors[3]=ID+1

        # Southeast
        if (self.tot_size-ID< self.num_cols+1):
            Neighbors[4]= -1
        elif ((ID+1)%self.num_cols==0):
            Neighbors[4]= -1
        else:
            Neighbors[4]= ID+1+self.num_cols

        # South
        if (self.tot_size-ID < self.num_cols +1):
            Neighbors[5]=-1
        else:
            Neighbors[5]= ID+ self.num_cols

        # Southwest
        if (self.tot_size-ID < self.num_cols+1):
            Neighbors[6]= -1
        elif (ID% self.num_cols==0):
            Neighbors[6]=-1
        else:
            Neighbors[6]= ID+self.num_cols-1

        # West
        if (ID+1)%self.num_cols==0:
            Neighbors[7]= -1
        else:
            Neighbors[7]= ID-1

        return Neighbors

    def _flow_dir(self, ID_from: int, ID_to: int) -> int:
        '''
        This function sets the flow direction according to the ID.
        flow direction is from current cell to cell which caused it to flood, clockwise from East (1,2,4,8,16,32,64,128)
        If two cells are not neighbors or if a neighbor is off the grid/ has no_data, direction set to 0.
        '''

        direction= np.zeros(8)
        if (ID_to== ID_from+1):
            #flow to east
            direction= 1
        elif (ID_to== ID_from+self.num_cols+1):
            #flow to southeast
            direction=2
        elif (ID_to== ID_from+self.num_cols):
            #flow to south
            direction=4
        elif (ID_to== ID_from +self.num_cols-1):
            #flow to southwest
            direction=8
        elif (ID_to== ID_from -1):
            #flow to west
            direction=16
        elif (ID_to== ID_from -self.num_cols-1):
            #flow to northwest
            direction=32
        elif (ID_to== ID_from -self.num_cols):
            #flow to north
            direction=64
        elif (ID_to== ID_from -self.num_cols+1):
            #flow to northeast
            direction= 128
        else:
            #cells are not neighbors
            direction= 0

        return direction

    def _repr_neighbors(self, ID, neighbors):
        '''
        This function is the representation of the neighbors of one cell according to the Neighbors vector
        _______
        |1|2|3|
        |4|5|6|
        |7|8|9|
        _______
        '''
        mat= np.zeros((3,3), dtype=np.int16)
        mat[1,1]= ID
        mat[0,:]=neighbors[:3]
        mat[1,2]=neighbors[3]
        mat[2,:]= neighbors[6:3:-1]
        mat[1,0]=neighbors[7]

        return mat

    def isLocalMinimum(self, ID: int) -> bool:
        minimum= True
        m,n= ID//self.num_cols, ID% self.num_rows
        if self.geo_info['flooded'][m,n]==2: #current cell is on confirmed path to outlet
            minimum= False
        else:
            neighbors= self._get_neighbors(ID)
            while i<8 and minimum==True:
                if neighbors[i]>-1:
                    neighbor_m, neighnor_n= neighbors[i]//self.num_cols, neighbors % self.num_rows
                    if self.geo_info['dem'][neighbor_m,neighbor_n]<=self.geo_info['dem'][m,n]:
                        minimum= False
                i+=1

        return minimum

    def _id_transform(self, ID):
        '''Simply transform ID to rows and columns value'''
        return int(ID//self.num_cols), int(ID % self.num_rows)

    def _read_dem_layer(self, layer):
        ''' Input layer: QGIS loaded project layer'''
        gd= gdal.Open(layer.source())
        array= gd.ReadAsArray()
        #get layer information
        proj= gd.GetProjection()
        geo_trans= gd.GetGeoTransform()
        rows, cols= array.shape
        bands= 1

        self.metadata= {'projection': proj,
                        'geotransform': geo_trans,
                        'rows': rows,
                        'cols': cols,
                        'bands': bands,
                        }

        return array

    def _read_dem_tif(self, dem_map):
        '''return dem data, bounds, crs'''
        raster= rasterio.open(dem_map)
        data= raster.read(1) #in case raster data only has 1 band
        self.map_info= {'crs': raster.crs,
                        'bounds': raster.bounds,
                        'height': raster.height,
                        'width': raster.width,
                        'count': raster.count,
                        'transform': raster.transform,}

        return data

    def _write_dem_gdal(self, src=None, dst=None):
        '''
        src: dem array
        dst: str, point to where you want to store the file, if not specified, create a temporary repo instead.
        '''
        if not os.path.exists('temp') and dst is None:
            os.mkdir('temp')
            dst= os.path.join('temp','test.tif')
        if src is None:
            src= self.geo_info['dem']
        driver= gdal.GetDriverByName("GTiff")
        rows= self.metadata.get('rows')
        cols= self.metadata.get('cols')
        geo_trans= self.metadata.get('geotransform')
        bands= self.metadata.get('bands')
        crs= self.metadata.get('projection')
        outdata = driver.Create(dst, rows, cols, bands, gdal.GDT_Float32)
        outdata.SetGeoTransform(geo_trans)
        outdata.SetProjection(crs)
        outdata.GetRasterBand(1).WriteArray(src)
        outdata.FlushCache()
        outdata = None

        return None

    def _write_dem(self, src= None, dst= None, height= None, width= None,
                     crs=None, transform=None, map=None):
        if self.map_info is not None:
            map= self.map_info
        src= self.geo_info.get('dem', None)
        height= map.get('height', None)
        width= map.get('width', None)
        crs= map.get('crs', None)
        transform= map.get('transform', None)
        count= map.get('count', None)
        dtype= map.get('dtype', float)
        if src is None or crs is None or height is None or dst is None:
            raise ValueError("Input argument incomplete!")
        with rasterio.open(dst, 'w',
                          driver='GTiff',
                          height=height,
                          width=width,
                          count=1,
                          dtype=src.dtype,
                          crs=crs,
                          transform=transform,
                          ) as f:
            f.write(src, 1)

class Point(object):
    def __init__(self, ID=0, elev=0, name=None):
        self.id= ID
        self.elev= elev
        self.name= name

    def __repr__(self):
        return 'name:    %s; ID:    %d; elevation:    %.2f'%(self.name, self.ID, self.elevation)

class Queue(object):
    def __init__(self, name):
        self.id= []
        self.elev= []

    def __repr__(self):
        return self.name

    def push(self, pnt):
        #pnt: Point object
        self.id.append(pnt.id)
        self.elev.append(pnt.elev)

    def top(self):
        # retrieve the smallest num
        id_arr= np.array(self.id)
        elev_arr= np.array(elev_arr)
        self.del_ind= np.argmin(elev_arr)
        cur_id= id_arr[ind]
        pnt= Point(cur_id, elev_arr.min())

        return pnt, ind

    def pop(self):
        self.id.pop(self.del_ind)
        self.elev.pop(self.del_ind)
        del self.del_ind

    def empty(self):
        return bool(self.id)
