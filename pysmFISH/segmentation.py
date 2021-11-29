import pickle
import itertools
import numpy as np





def load_stitched_segmented_data(fov,segmentation_output_path):
    segmented_region_img = pickle.load(open(segmentation_output_path / ('registered_objs_dict_fov_' + str(fov) + '.pkl'), 'rb'))

    obj_segmented_img = {el:coords_dict['stitched_coords'] for (el,coords_dict) in segmented_region_img.items()}
    return obj_segmented_img

def objects_overlapping_region(obj_segmented_img,chunk_coords):
    objs_list=list()
    for lab in obj_segmented_img.keys():
        if ((obj_segmented_img[lab][:,0] >= chunk_coords[0]).any() and
            (obj_segmented_img[lab][:,0] <= chunk_coords[1]).any() and
            (obj_segmented_img[lab][:,1] >= chunk_coords[2]).any() and
            (obj_segmented_img[lab][:,1] <= chunk_coords[3]).any()):
                objs_list.append(lab)
    return objs_list

def overlapping_objs(obj_segmented_img1,obj_segmented_img2,chunk_coords,min_overlapping_pixels_segmentation = 20):
    objs_list_img1 = objects_overlapping_region(obj_segmented_img1,chunk_coords)
    objs_list_img2 = objects_overlapping_region(obj_segmented_img2,chunk_coords)
    
    all_combinations = list(itertools.product(objs_list_img1,objs_list_img2))
    
    intersection = []
    overlapping_sizes = {}
    for couple in all_combinations:
        obj1_set = set(map(tuple, obj_segmented_img1[couple[0]].astype(int)))
        obj2_set = set(map(tuple, obj_segmented_img2[couple[1]].astype(int)))
        result=not obj1_set.isdisjoint(obj2_set)
        if result:
            # Check minimum number of pixels to correct for potential error caused by the rounding
            number_overlapping_pixels = len(obj1_set.intersection(obj2_set))
            if number_overlapping_pixels >= min_overlapping_pixels_segmentation:
                intersection.append(result)
                overlapping_sizes[couple] = number_overlapping_pixels
            else:
                intersection.append(not result)
        else:
            intersection.append(result)
    if (np.array(intersection)).any():
        all_combinations = np.array(all_combinations)
        intersecting_cpls = all_combinations[intersection]
    else:
        intersecting_cpls = []
    return intersecting_cpls, overlapping_sizes



def correct_multiple_overlapping_objects(intersecting_cpls,overlapping_sizes):

    intersecting_cpls_tpl_trimmed = []
    intersecting_cpls_tpl = [tuple(x) for x in intersecting_cpls]
    
    intersecting_cpls_removed_overlapping = intersecting_cpls_tpl.copy()
    tracker = intersecting_cpls_tpl.copy()

    combined_intersecting_cpls_tpl = []
    tracker = []
    reducer = intersecting_cpls_tpl.copy()
    for tpl in intersecting_cpls_tpl:
        if tpl not in tracker:
            combined = [tpl]
            tracker.append(tpl)
            reducer.remove(tpl)
            runner = reducer
            for cpl in runner:
                if set(tpl).intersection(cpl):
                    tracker.append(cpl)
                    combined.append(cpl)
                    reducer.remove(cpl)
            combined_intersecting_cpls_tpl.append(combined)

    #combined_intersecting_cpls_tpl = set(combined_intersecting_cpls_tpl)
    for grp in combined_intersecting_cpls_tpl:
        if len(grp) >1 :
            # Select the object with largest overlap (maybe not best criteria)
            selected = { your_key: overlapping_sizes[your_key] for your_key in grp }
            counts = list(selected.values())
            keys = list(selected.keys())
            selected_cpl = keys[counts.index(max(counts))]
            intersecting_cpls_tpl_trimmed.append(selected_cpl)
        else:
            intersecting_cpls_tpl_trimmed.append(grp)
    
    return intersecting_cpls_removed_overlapping,overlapping_sizes



def merge_objects(obj_segmented_img1,obj_segmented_img2,intersecting_cpls):
    merged_objs = {}
    for cpl in intersecting_cpls:
        obj1 = obj_segmented_img1[cpl[0]].astype(int)
        obj2 = obj_segmented_img2[cpl[1]].astype(int)
        merged_obj = np.vstack((obj1,obj2))
        merged_obj = np.unique(merged_obj,axis=0)
        merged_objs[cpl[0]] = merged_obj
    return merged_objs



def adjust_overlapping_objects(obj_segmented_img1,obj_segmented_img2,chunk_coords,segmentation_output_path,min_overlapping_pixels_segmentation = 20):
    
    intersecting_cpls, overlapping_sizes = overlapping_objs(obj_segmented_img1,obj_segmented_img2,chunk_coords,min_overlapping_pixels_segmentation)
    
    # Correct for multiple overlapping large objects
    intersecting_cpls,overlapping_sizes = correct_multiple_overlapping_objects(intersecting_cpls,overlapping_sizes)
    
    objs = merge_objects(obj_segmented_img1,obj_segmented_img2,intersecting_cpls)
    
    #for cpl in intersecting_cpls:
     #   del obj_segmented_img1[cpl[0]]
     #   del obj_segmented_img2[cpl[1]]
     #   obj_segmented_img1[cpl[0]] = objs[cpl[0]]
    
    return intersecting_cpls, objs



def remove_overlapping_obj(fov,
                           nuclei_org_tiles,
                          segmentation_output_path,
                          min_overlapping_pixels_segmentation = 20):

    overlapping_dict = nuclei_org_tiles.overlapping_regions[fov]
    all_obj_segmented = {}
    output_cpl = {}
    corner_cpl = []
    
    
    # Correct overlapping objects at the intersection between the images
    for cpl, chunk_coords in overlapping_dict.items():
        obj_segmented_img1 = load_stitched_segmented_data(cpl[0],segmentation_output_path)
        obj_segmented_img2 = load_stitched_segmented_data(cpl[1],segmentation_output_path)
        intersecting_cpls, objs = adjust_overlapping_objects(obj_segmented_img1,obj_segmented_img2,chunk_coords,segmentation_output_path,min_overlapping_pixels_segmentation)
        all_obj_segmented[cpl[0]] = obj_segmented_img1
        all_obj_segmented[cpl[1]] = obj_segmented_img2
        output_cpl[cpl] = (intersecting_cpls, objs)
        corner_cpl.append(cpl[1])

        
    # Check if processing fov has one single neighbour
    if len(overlapping_dict.keys()) == 2:
        # Correct overlapping objects at the corner
        intersecting_cpls, objs = adjust_overlapping_objects(all_obj_segmented[corner_cpl[0]],all_obj_segmented[corner_cpl[1]],chunk_coords,segmentation_output_path,min_overlapping_pixels_segmentation)
        all_obj_segmented[corner_cpl[0]] = obj_segmented_img1
        all_obj_segmented[corner_cpl[1]] = obj_segmented_img2
        output_cpl[tuple(corner_cpl)] = (intersecting_cpls, objs)
    
    
    # Format the output
    all_cpls = []
    all_intersecting_step = []
    all_obj_step = {}
    for cpl, (intersecting_cpls, objs) in output_cpl.items():
        all_intersecting_step.append(intersecting_cpls)
        all_obj_step.update(objs)
        all_cpls.append(cpl)

    unraveled_intersecting_step = [val for big_list in all_intersecting_step for el in big_list for val in el ]
    
    # To use in case to sort by fov for parallelisation
    #all_couples = set([el for cpl in all_cpls for el in cpl])
    #for el in all_couples:
     #   intersecting_step_dict[el] = [val for val in unraveled_intersecting_step if str(el) in val]

    
    return unraveled_intersecting_step,all_obj_step