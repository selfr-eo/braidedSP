# Typing imports
from dataclasses import dataclass

# general imports
from datetime import datetime

# functional imports
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, MultiPoint, LineString




@dataclass
class Centerline:

    river_name: str
    date: datetime
    gdf: gpd.GeoDataFrame

    def __post_init__(self):

        self.crs = self.gdf.crs
    
    @classmethod
    def from_file(cls, river_name, date, path):

        gdf = gpd.read_file(path)

        return cls(river_name, date, gdf)



    def trim_to_river_bounds(self, river_bounds):

        """Depending on the watermask and starting point of the river system, there can be circular loops generated at the start of the system.
        By trimming the centerline to a start line (river cross section) we can ensure that the starting reach(es) is(are) clean.

        Parameters
        ----------
        cl : _type_
            _description_
        starting_line : _type_
            _description_
        """
        # make sure we are working in the crs of the centerlines
        river_bounds = river_bounds.to_crs(self.gdf.crs)

        # cut the geodataframe by the centerline
        cl_new = self.gdf.copy()
        cl_new = gpd.clip(cl_new, river_bounds.unary_union)

        return cl_new



    def merge_short_centerlines(self, connection_threshold=50, similarity_threshold=0.7, shortest_centerline=2000):

        # 1. Ensure centerlines are in utm crs
        utm_crs = self.gdf.estimate_utm_crs()
        cl_utm = self.gdf.to_crs(utm_crs)

        # Calculate the length of each branch
        cl_utm['length'] = cl_utm.length

        # 2. Iteratively merge branches
        while True:

            # extract all branches as a geodataframe that are less than 2 km
            short_branches = cl_utm.loc[cl_utm['length'] < shortest_centerline]

            # Stopping condition is when we have no more branches remaining
            if short_branches.empty:
                break

            # 
            for branch_id, branch_geom in short_branches[['branch_id', 'geometry']].values:
                branch_endpoints = _get_endpoints_and_vectors(branch_geom)
                closest_candidates = []

                # Find connecting branches within the threshold distance
                for candidate_id, candidate_geom in cl_utm[['branch_id', 'geometry']].values:
                    if branch_id == candidate_id:
                        continue

                    candidate_endpoints = _get_endpoints_and_vectors(candidate_geom)

                    for branch_endpoint, branch_vector in branch_endpoints:
                        for candidate_endpoint, candidate_vector in candidate_endpoints:
                            distance = Point(branch_endpoint).distance(Point(candidate_endpoint))
                            if distance <= connection_threshold:
                                similarity = cosine_similarity(
                                    branch_vector.reshape(1, -1), 
                                    candidate_vector.reshape(1, -1)
                                )[0, 0]
                                closest_candidates.append({
                                    "branch_id": candidate_id,
                                    "distance": distance,
                                    "similarity": similarity,
                                    "branch_endpoint": branch_endpoint,
                                    "candidate_endpoint": candidate_endpoint
                                })

                # Select the best candidate based on similarity
                if closest_candidates:
                    best_candidate = max(
                        closest_candidates,
                        key=lambda x: (x['similarity'], -x['distance'])  # Prioritize similarity, then proximity
                    )
                    if best_candidate['similarity'] > similarity_threshold:
                        # Merge branches
                        candidate_id = best_candidate["branch_id"]
                        candidate_geom = cl_utm.loc[cl_utm['branch_id'] == candidate_id, 'geometry'].iloc[0]

                        # Combine coordinates
                        all_coords = []
                        if best_candidate['branch_endpoint'] == tuple(branch_geom.coords[0]):
                            all_coords.extend(list(branch_geom.coords[::-1]))
                        else:
                            all_coords.extend(list(branch_geom.coords))

                        if best_candidate['candidate_endpoint'] == tuple(candidate_geom.coords[0]):
                            all_coords.extend(list(candidate_geom.coords))
                        else:
                            all_coords.extend(list(candidate_geom.coords[::-1]))

                        # Update GeoDataFrame
                        merged_geom = LineString(all_coords)
                        cl_utm = cl_utm[cl_utm['branch_id'] != branch_id]
                        cl_utm = cl_utm[cl_utm['branch_id'] != candidate_id]

                        new_row = pd.DataFrame({
                            "branch_id": [candidate_id],  # Retain candidate ID (the larger branch)
                            "geometry": [merged_geom],
                            "length": [merged_geom.length]
                        })
                        cl_utm = pd.concat([cl_utm, gpd.GeoDataFrame(new_row, crs=cl_utm.crs)])
                        break  # Restart the loop after updating

            else:
                # If no merges occurred, exit the loop
                break

        cl_new  = cl_utm.sort_values(by=['branch_id'])
        cl_new = cl_new.to_crs(self.crs)
        cl_new = cl_new.reset_index(drop=True)
        cl_new['branch_id_old'] = cl_new['branch_id']
        cl_new['branch_id'] = cl_new.index + 1


        # Return the updated GeoDataFrame
        return cl_new



    def join_cl_at_joints(self, starting_line, main_centerline, search_dist=20):

        """Reaches are seperated due to the processing of pixels and selecting their center coordinates. we join them here and also keep track of which reach each reach flows into"""

        # fist convert centerline to local coordinates
        cl = self.gdf.copy()
        local_crs = cl.estimate_utm_crs()
        cl = cl.to_crs(local_crs)
        starting_line = starting_line.to_crs(local_crs)
        main_centerline = main_centerline.to_crs(local_crs)
        # main_centerline = linemerge(main_centerline.unary_union)
        main_centerline = main_centerline.geometry[0]


        # see if we need to reverse coordinates in main centerline
        main_centerline_endpoints = main_centerline.coords[0], main_centerline.coords[-1]
        dist = [Point(coord).distance(starting_line.unary_union) for coord in main_centerline_endpoints]
        if dist[0] > dist[-1]: # reverse
            reversed_coords = list(main_centerline.coords)
            reversed_coords.reverse()
            main_centerline = LineString(reversed_coords)

        
        # Cycle through centerlines, and document endpoints
        branch_inlet = list()
        branch_outlet = list()
        branch_proj_dist = list()
        for i in cl.index:

            branch_id = cl.loc[i, 'branch_id']
            line = cl.loc[i, 'geometry']

            # extract coords, and order from first to last according to main center line
            dist   = [main_centerline.project( Point(coord), normalized=True ) for coord in [line.coords[0], line.coords[-1]]]

            # see if we need to reverse the line string
            if dist[0] > dist[-1]: # reverse
                reversed_coords = list(line.coords)
                reversed_coords.reverse()
                line = LineString(reversed_coords)
                cl.loc[i, 'geometry'] = line

            # save endpoints and assign the reach start dist and the reach dist
            inlet = Point(line.coords[0])
            outlet = Point(line.coords[-1]) # order the points
            proj_dist = np.min(dist)

            # add values to list
            branch_inlet.append(inlet)
            branch_outlet.append(outlet)
            branch_proj_dist.append(proj_dist)

        # add to dataframe, create new index and sort by distance
        inlet_coords = dict(zip(cl.branch_id.values, branch_inlet))
        outlet_coords = dict(zip(cl.branch_id.values, branch_outlet))
        cl['branch_proj_dist'] = branch_proj_dist
        cl = cl.sort_values(by='branch_proj_dist').reset_index(drop=True)

        # now cycle through dataframe and assign downstream ids, only assigning to upstream points
        cl['downstream_ids'] = None
        for i in cl.index:

            branch_id = cl.loc[i, 'branch_id']

            # Check what inlet coords are within the search distance if the branches outlet
            downstream_branches = list()
            for key in inlet_coords.keys():
                if key != branch_id:
                    if outlet_coords[branch_id].distance(inlet_coords[key]) < search_dist:
                        downstream_branches.append(key)

            # Check also if there are neighboring outlet branches that should be joined
            shared_outlet_branches = list()
            for key in inlet_coords.keys():
                if key != branch_id:
                    if outlet_coords[branch_id].distance(outlet_coords[key]) < search_dist:
                        shared_outlet_branches.append(key)

            # add the list of downstram brach ids to the dataframe
            cl.loc[i, 'downstream_ids'] = ', '.join([str(key) for key in downstream_branches])

            # Find the new common jon for the outlet and downstream inlets
            shared_point = MultiPoint([inlet_coords[key] for key in downstream_branches]+ [outlet_coords[key] for key in shared_outlet_branches] + [outlet_coords[branch_id]]).centroid
            outlet_coords[branch_id] = shared_point
            for key in downstream_branches:
                inlet_coords[key] = shared_point
            for key in shared_outlet_branches:
                outlet_coords[key] = shared_point

        # Finally, update endpoints to the revised editions
        # loop through each branch, and update the endpoints with whats in the updated dictionaries
        for i in cl.index:

            # extract points in linestring for branch
            branch_id = cl.loc[i, 'branch_id']
            coords = [Point(coord) for coord in cl.loc[i, 'geometry'].coords]
            coords[0] = inlet_coords[branch_id]
            coords[-1] = outlet_coords[branch_id]

            # update in place
            cl.loc[i, 'geometry'] = LineString(coords)

        return cl
    

    def export(self, path):
        self.cl.to_file(path)




def _calculate_vector(branch_geom, endpoint_index):
    # Calculate the direction vector for a LineString at the given endpoint.
    coords = list(branch_geom.coords)
    if endpoint_index == 0:  # Start of the branch
        vec = np.array(coords[1]) - np.array(coords[0])
    else:  # End of the branch
        vec = np.array(coords[-1]) - np.array(coords[-2])
    return vec / np.linalg.norm(vec)  # Normalize


def _get_endpoints_and_vectors(linestring):
    """Return the endpoints and their respective vectors."""
    return [
        (tuple(linestring.coords[0]), _calculate_vector(linestring, 0)),  # Start point
        (tuple(linestring.coords[-1]), _calculate_vector(linestring, -1)),  # End point
    ]
