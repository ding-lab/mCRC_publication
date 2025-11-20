import Morph
import numpy
import csv
import matplotlib.pyplot
import skimage
import pandas as pd
import os
import time
import pickle
import scipy
from Morph.features import Distance

# Configuration
# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_FILE_PATH = os.path.join(SCRIPT_DIR, 'mCRC_Morph_Mask_tracking.csv')
OUTDIR = os.path.join(SCRIPT_DIR, "Out")
MORPH_ANNOTATION_PATH = os.path.join(OUTDIR, "xenium_morph_annotations.pkl")

# Load data from CSV file instead of Google Sheets
print(f"Loading data from CSV file: {CSV_FILE_PATH}")
mask_df = pd.read_csv(CSV_FILE_PATH)
print(f"Loaded {len(mask_df)} rows from CSV")

# Helper function to safely get string values from CSV rows (handles NaN)
def safe_str(value):
    """Convert value to string, handling NaN/None values."""
    if pd.isna(value) or value == 'NULL' or value == '':
        return None
    return str(value).strip()

# Helper function to safely split comma-separated values
def safe_split(value, default=None):
    """Split comma-separated string, handling NaN/None values."""
    if pd.isna(value) or value == 'NULL' or value == '':
        return default if default is not None else []
    return str(value).strip().split(',')

# Helper function to create a clean set from comma-separated cluster IDs
def clean_cluster_set(value):
    """Create a set of cluster IDs from comma-separated string, removing empty values."""
    if pd.isna(value) or value == 'NULL' or value == '':
        return set()
    parts = str(value).strip().split(',')
    # Remove empty strings and strip whitespace
    cluster_set = {part.strip() for part in parts if part.strip()}
    return cluster_set

# Helper function to check if any cluster IDs in G exist in the data
def has_matching_clusters(data, G):
    """Check if any cluster IDs in G exist in the mapped data."""
    if 'g' not in data or len(data['g']) == 0:
        return False
    unique_clusters = set(data['g'])
    return bool(G & unique_clusters)

# 0. Global Functions

def clusters(file, cluster_name):
    data = {}
    with open(file, 'rt') as f:
        dict_reader = csv.DictReader(f)
        for row in dict_reader:
            data[row['cell_id']] = row[cluster_name]
    return data

def map(clusters, cells):
    g = []
    x = []
    y = []
    v = []
    for c, a, b in zip(cells['g'], cells['x'], cells['y']):
        if c in clusters:
            g.append(clusters[c])
            x.append(a)
            y.append(b)
            v.append(1)
    return {'g': g, 'x': numpy.array(x), 'y': numpy.array(y), 'v': numpy.array(v)}

def annotation(file, shape, d):
    points = []
    with open(file, 'rt') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            i += 1
            if i < 4:
                continue
            points.append((float(row[1]) / d, float(row[2]) / d))
    image = skimage.draw.polygon2mask(shape, points)
    return image

def hole_annotation(file, shape, d):
    points = []
    with open(file, 'rt') as f:
        reader = csv.reader(f)
        i = 0
        selection = 'Selection 1'
        image = numpy.zeros(shape)
        for row in reader:
            i += 1
            if i < 4:
                continue
            if row[0] != selection: 
                selection = row[0]
                image += skimage.draw.polygon2mask(shape, points)
                points = []
            points.append((float(row[1]) / d, float(row[2]) / d))
        # Add the mask for the final selection.
        if points:
            image += skimage.draw.polygon2mask(shape, points)
    return image

# Load existing annotations if available
if os.path.exists(MORPH_ANNOTATION_PATH):
    print(f"Loading existing annotations from {MORPH_ANNOTATION_PATH}")
    with open(MORPH_ANNOTATION_PATH, "rb") as f:
        xenium_annotations = pickle.load(f)
    print(f"Loaded annotations for {len(xenium_annotations)} sections")
else:
    print("No existing annotations found. Starting fresh.")
    xenium_annotations = {}

# 1. Read Data
print("\n=== Step 1: Reading Data ===")
xenium_annotations = {}

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    banksy_path = row['Banksy_annotation']  
    xenium_path = row['Xenium_output']      

    # Read cluster data
    print(f"Processing {section_id} - Clusters: {banksy_path}")
    cluster_data = clusters(banksy_path, 'banksy_leiden_0.7')

    valid_cell_ids = set(cluster_data.keys())

    # Read cell data
    print(f"Processing {section_id} - Cells: {xenium_path}")
    cell_data = Morph.readers.cells(xenium_path + '/cells.csv.gz') 

    filtered_cell_data = {'g': [], 'x': [], 'y': [], 'v': []}
    for g, x, y in zip(cell_data['g'], cell_data['x'], cell_data['y']):
        if g in valid_cell_ids:
            filtered_cell_data['g'].append(g)
            filtered_cell_data['x'].append(x)
            filtered_cell_data['y'].append(y)
            filtered_cell_data['v'].append(1) 
    filtered_cell_data['x'] = numpy.array(filtered_cell_data['x'])
    filtered_cell_data['y'] = numpy.array(filtered_cell_data['y'])
            
    xenium_annotations[section_id] = {}
    
    # Apply map function
    if cluster_data and filtered_cell_data is not None:
        xenium_annotations[section_id]['cell'] = filtered_cell_data
        mapped_data = map(cluster_data, filtered_cell_data)
        xenium_annotations[section_id]['data'] = mapped_data
        print(f"Completed {section_id}: {len(mapped_data['g'])} cells mapped")
    else:
        print(f"Skipping {section_id} due to missing data")
        xenium_annotations[section_id] = None

# 2. NAT mask
print("\n=== Step 2: NAT mask ===")

## A. Mucosa mask
def mucosa_area_filter(image, lambda_area_closing, hole_path):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing).astype(int)
    hole_mask = hole_annotation(hole_path, image.shape, d)
    output_image = ((image + hole_mask) > 0).astype(int)    
    return output_image

d = 10
G = {'4', '7', '14', '17'}
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing Colon mucosa mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)
    
    colon_mucosa_coords = safe_str(row['Colon_mucosa_coordinates'])
    if colon_mucosa_coords is not None:
        hole_path = safe_str(row['CM_holes_coordinates'])
        if hole_path is None:
            raise ValueError(f"CM_holes_coordinates is required when Colon_mucosa_coordinates is provided for {section_id}")
        
        image = Morph.backbone(xenium_annotations[section_id]['data'], 
                               ['xenium', d], 
                               ['total', G], 
                               ['maximum'], 
                               ['naive'], 
                               ['binary', tau],  
                               ['custom', mucosa_area_filter, 100000, hole_path], 
                               ['blob', S])
        file_path = colon_mucosa_coords
        mask = annotation(file_path, image.shape, d)
        cm_mask = numpy.multiply(image > 0, mask).astype(int)
        
        xenium_annotations[section_id]['cm_mask'] = cm_mask
        
        mapper = Morph.modules.Mapper()
        cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

        output_file = os.path.join(section_dir, f"{section_id}_cm_mask.csv")
        Morph.writers.xenium(output_file, cm_mask, cells)
        print(f"write {output_file}")
    else:
        cm_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                               int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))

        xenium_annotations[section_id]['cm_mask'] = cm_mask
        
        mapper = Morph.modules.Mapper()
        cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

        output_file = os.path.join(section_dir, f"{section_id}_cm_mask.csv")
        Morph.writers.xenium(output_file, cm_mask, cells)
        print(f"write {output_file}")

## B. Colon lymphoid follicle mask
### Functions
def cl_filter(image, mask, element=None):
    image = numpy.multiply(image, mask)
    return Morph.operators.opening(Morph.operators.closing(image, element), element)

### Post-processing
d = 10
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing Colon lymphoid mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    colon_lymphoid = safe_str(row['Colon_lymphoid'])
    if colon_lymphoid is None:
        cl_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                              int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
        
        xenium_annotations[section_id]['cl_mask'] = cl_mask
        
        mapper = Morph.modules.Mapper()
        cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

    else:
        G = clean_cluster_set(colon_lymphoid)
        if not G:
            print(f"Warning: Colon_lymphoid is empty for {section_id}, creating empty mask")
            cl_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                                  int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
            xenium_annotations[section_id]['cl_mask'] = cl_mask
            mapper = Morph.modules.Mapper()
            cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
            continue
        
        # Check if any clusters in G exist in the data
        if not has_matching_clusters(xenium_annotations[section_id]['data'], G):
            print(f"Warning: No matching Colon_lymphoid clusters in data for {section_id}, creating empty mask")
            cl_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                                  int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
            xenium_annotations[section_id]['cl_mask'] = cl_mask
            mapper = Morph.modules.Mapper()
            cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
            continue
    
        image = Morph.backbone(xenium_annotations[section_id]['data'], 
                               ['xenium', d], 
                               ['total', G], 
                               ['maximum'],  
                               ['custom', cl_filter, 1 - xenium_annotations[section_id]['cm_mask'], S], 
                               ['binary', tau], 
                               ['naive'], 
                               ['blob', S])

        cl_to_keep = safe_split(row['CL_to_keep'])
        if cl_to_keep:
            keep_values = [int(i) for i in cl_to_keep]
            image[~numpy.isin(image, keep_values)] = 0

        mapper = Morph.modules.Mapper()
        cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

        cl_label = image
        xenium_annotations[section_id]['cl_label'] = cl_label

        output_file = os.path.join(section_dir, f"{section_id}_cl_label.csv")
        Morph.writers.xenium(output_file, cl_label, cells)
        print(f"write {output_file}")

    
        cl_mask = (image > 0).astype(int)
        xenium_annotations[section_id]['cl_mask'] = cl_mask

## C. Liver parenchyma mask
### Functions
def lp_area_filter(image, lambda_area_closing, lambda_area_opening):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    return skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)

def lp_hole_area_filter(image, lambda_area_closing, lambda_area_opening, hole_path):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    image = skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)
    hole_mask = hole_annotation(hole_path, image.shape, d)
    output_image = ((image + hole_mask) > 0).astype(int)
    return output_image

### Post-processing
d = 10
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing Liver parenchyma mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    liver_parenchyma = safe_str(row['Liver_parenchyma'])
    if liver_parenchyma is None:
        lp_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                              int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
        
        xenium_annotations[section_id]['lp_mask'] = lp_mask
        
        mapper = Morph.modules.Mapper()
        cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

    else:
        nat_holes = safe_str(row['NAT_holes_coordinates'])
        if nat_holes is not None:
            hole_path = nat_holes
            
            if section_id != 'HT413C1-Th1K2A1U2':
                
                G = clean_cluster_set(liver_parenchyma)
                if not G:
                    print(f"Warning: Liver_parenchyma is empty for {section_id}, creating empty mask")
                    lp_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                                          int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
                    xenium_annotations[section_id]['lp_mask'] = lp_mask
                    mapper = Morph.modules.Mapper()
                    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
                    continue
                
                # Check if any clusters in G exist in the data
                if not has_matching_clusters(xenium_annotations[section_id]['data'], G):
                    print(f"Warning: No matching Liver_parenchyma clusters in data for {section_id}, creating empty mask")
                    lp_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                                          int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
                    xenium_annotations[section_id]['lp_mask'] = lp_mask
                    mapper = Morph.modules.Mapper()
                    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
                    continue
        
                image = Morph.backbone(xenium_annotations[section_id]['data'], 
                                       ['xenium', d], 
                                       ['total', G], 
                                       ['maximum'],  
                                       ['naive'], 
                                       ['binary', tau], 
                                       ['custom', lp_hole_area_filter, 100000, 2, hole_path],  
                                       ['blob', S])

                lp_to_keep = safe_split(row['LP_to_keep'])
                if lp_to_keep:
                    keep_values = [int(i) for i in lp_to_keep]
                    image[~numpy.isin(image, keep_values)] = 0
    
                mapper = Morph.modules.Mapper()
                cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
    
                lp_label = image
                xenium_annotations[section_id]['lp_label'] = lp_label
    
                output_file = os.path.join(section_dir, f"{section_id}_lp_label.csv")
                Morph.writers.xenium(output_file, lp_label, cells)
                print(f"write {output_file}")
    
        
                lp_mask = (image > 0).astype(int)
                xenium_annotations[section_id]['lp_mask'] = lp_mask
            else:
                G = clean_cluster_set(liver_parenchyma)
                if not G:
                    print(f"Warning: Liver_parenchyma is empty for {section_id}, creating empty mask")
                    lp_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                                          int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
                    xenium_annotations[section_id]['lp_mask'] = lp_mask
                    mapper = Morph.modules.Mapper()
                    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
                    continue
                
                # Check if any clusters in G exist in the data
                if not has_matching_clusters(xenium_annotations[section_id]['data'], G):
                    print(f"Warning: No matching Liver_parenchyma clusters in data for {section_id}, creating empty mask")
                    lp_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                                          int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
                    xenium_annotations[section_id]['lp_mask'] = lp_mask
                    mapper = Morph.modules.Mapper()
                    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
                    continue
                
                image = Morph.backbone(xenium_annotations[section_id]['data'], 
                                       ['xenium', d], 
                                       ['total', G], 
                                       ['maximum'],  
                                       ['naive'], 
                                       ['binary', tau], 
                                       ['custom', lp_hole_area_filter, 10, 2, hole_path],  
                                       ['blob', S])
                
                lp_to_keep = safe_split(row['LP_to_keep'])
                if lp_to_keep:
                    keep_values = [int(i) for i in lp_to_keep]
                    image[~numpy.isin(image, keep_values)] = 0
    
                mapper = Morph.modules.Mapper()
                cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
    
                lp_label = image
                xenium_annotations[section_id]['lp_label'] = lp_label
    
                output_file = os.path.join(section_dir, f"{section_id}_lp_label.csv")
                Morph.writers.xenium(output_file, lp_label, cells)
                print(f"write {output_file}")
    
        
                lp_mask = (image > 0).astype(int)
                xenium_annotations[section_id]['lp_mask'] = lp_mask

## D. Lung & Breast parenchyma mask
def lup_hole_area_filter(image, lambda_area_closing, hole_path):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    hole_mask = hole_annotation(hole_path, image.shape, d)
    output_image = ((image + hole_mask) > 0).astype(int)
    return output_image

### Processing
d = 10
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing Lung & Breast parenchyma mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    lubr_parenchyma = safe_str(row['LuBr_parenchyma'])
    if lubr_parenchyma is None:
        lubnat_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                              int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
        
        xenium_annotations[section_id]['lubnat_mask'] = lubnat_mask

    else:
        
        G = clean_cluster_set(lubr_parenchyma)
        if not G:
            print(f"Warning: LuBr_parenchyma is empty for {section_id}, creating empty mask")
            lubnat_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                              int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
            xenium_annotations[section_id]['lubnat_mask'] = lubnat_mask
            continue
        
        # Check if any clusters in G exist in the data
        if not has_matching_clusters(xenium_annotations[section_id]['data'], G):
            print(f"Warning: No matching LuBr_parenchyma clusters in data for {section_id}, creating empty mask")
            lubnat_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                              int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
            xenium_annotations[section_id]['lubnat_mask'] = lubnat_mask
            continue
        
        nat_holes = safe_str(row['NAT_holes_coordinates'])
        if nat_holes is not None:
            hole_path = nat_holes
            
            image = Morph.backbone(xenium_annotations[section_id]['data'], 
                                   ['xenium', d], 
                                   ['total', G], 
                                   ['maximum'],  
                                   ['naive'], 
                                   ['binary', tau], 
                                   ['custom', lup_hole_area_filter, 100, hole_path],  
                                   ['blob', S])
            
            lubr_to_keep = safe_split(row['LuBr_to_keep'])
            if lubr_to_keep:
                keep_values = [int(i) for i in lubr_to_keep]
                image[~numpy.isin(image, keep_values)] = 0

            mapper = Morph.modules.Mapper()
            cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

            lubnat_label = image
            xenium_annotations[section_id]['lubnat_label'] = lubnat_label

            output_file = os.path.join(section_dir, f"{section_id}_lubnat_label.csv")
            Morph.writers.xenium(output_file, lubnat_label, cells)
            print(f"write {output_file}")

    
            lubnat_mask = (image > 0).astype(int)
            xenium_annotations[section_id]['lubnat_mask'] = lubnat_mask

## E. Colon muscle mask
def colon_muscle_hole_area_filter(image, lambda_area_closing, lambda_area_opening, hole_path):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    image = skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int) 
    hole_mask = hole_annotation(hole_path, image.shape, d)
    output_image = ((image + hole_mask) > 0).astype(int)
    return output_image

d = 10
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing Colon muscle mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    colon_muscle = safe_str(row['Colon_Muscle'])
    if colon_muscle is None:
        c_mus_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                                     int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
        
        xenium_annotations[section_id]['c_mus_mask'] = c_mus_mask

    else:
        
        G = clean_cluster_set(colon_muscle)
        if not G:
            print(f"Warning: Colon_Muscle is empty for {section_id}, creating empty mask")
            c_mus_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                                     int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
            xenium_annotations[section_id]['c_mus_mask'] = c_mus_mask
            continue
        
        # Check if any clusters in G exist in the data
        if not has_matching_clusters(xenium_annotations[section_id]['data'], G):
            print(f"Warning: No matching Colon_Muscle clusters in data for {section_id}, creating empty mask")
            c_mus_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                                     int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
            xenium_annotations[section_id]['c_mus_mask'] = c_mus_mask
            continue
        
        c_mus_holes = safe_str(row['C_mus_holes_coordinates'])
        if c_mus_holes is not None:
            hole_path = c_mus_holes
            
            image = Morph.backbone(xenium_annotations[section_id]['data'], 
                                   ['xenium', d], 
                                   ['total', G], 
                                   ['maximum'],  
                                   ['naive'], 
                                   ['binary', tau], 
                                   ['custom', colon_muscle_hole_area_filter, 100, 50, hole_path],  
                                   ['blob', S])

            mapper = Morph.modules.Mapper()
            cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

            c_mus_label = image
            xenium_annotations[section_id]['c_mus_label'] = c_mus_label

            output_file = os.path.join(section_dir, f"{section_id}_c_mus_label.csv")
            Morph.writers.xenium(output_file, c_mus_label, cells)
            print(f"write {output_file}")

    
            c_mus_mask = (image > 0).astype(int)
            xenium_annotations[section_id]['c_mus_mask'] = c_mus_mask

# 3. Necrosis mask
print("\n=== Step 3: Necrosis mask ===")

## A. Functions
def n_filter(image, mask, element=None):
    image = numpy.multiply(image, mask)
    return Morph.operators.opening(Morph.operators.closing(image, element), element)

def n_area_filter(image, lambda_area_closing, lambda_area_opening):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    return skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)

## B. Processing
d = 10
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing Necrosis mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    necrosis = safe_str(row['Necrosis'])
    if necrosis is None:
        # Necrosis is NULL/NaN/empty - create empty mask
        n_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                              int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
        
        xenium_annotations[section_id]['n_mask'] = n_mask
        
        mapper = Morph.modules.Mapper()
        cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

        output_file = os.path.join(section_dir, f"{section_id}_n_mask.csv")
        Morph.writers.xenium(output_file, n_mask, cells)
        print(f"write {output_file}")

    else:
        G = clean_cluster_set(necrosis)
        if not G:
            print(f"Warning: Necrosis is empty for {section_id}, creating empty mask")
            n_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                              int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
            xenium_annotations[section_id]['n_mask'] = n_mask
            mapper = Morph.modules.Mapper()
            cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
            output_file = os.path.join(section_dir, f"{section_id}_n_mask.csv")
            Morph.writers.xenium(output_file, n_mask, cells)
            print(f"write {output_file}")
            continue
        
        # Check if any clusters in G exist in the data
        if not has_matching_clusters(xenium_annotations[section_id]['data'], G):
            print(f"Warning: No matching Necrosis clusters in data for {section_id}, creating empty mask")
            n_mask = numpy.zeros((int(xenium_annotations[section_id]['data']['x'].max() / d) + 1, 
                              int(xenium_annotations[section_id]['data']['y'].max() / d) + 1))
            xenium_annotations[section_id]['n_mask'] = n_mask
            mapper = Morph.modules.Mapper()
            cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)
            output_file = os.path.join(section_dir, f"{section_id}_n_mask.csv")
            Morph.writers.xenium(output_file, n_mask, cells)
            print(f"write {output_file}")
            continue

        n_size_filter = 0

        if section_id in ('CM732C1-AFR1U1', 'CM394C1-B1U1', 'CM626C1-A1U1', 'CM819C1-B3U1', 'CM426C1-A2U1'):
            n_size_filter = 30
    
        image = Morph.backbone(xenium_annotations[section_id]['data'], 
                               ['xenium', d], 
                               ['total', G], 
                               ['maximum'],  
                               ['custom', 
                                n_filter, 
                                1 - xenium_annotations[section_id]['cm_mask'] - xenium_annotations[section_id]['lp_mask'] - xenium_annotations[section_id]['lubnat_mask'], 
                                S], 
                               ['binary', tau], 
                               ['custom', n_area_filter, 10000, n_size_filter], 
                               ['blob', S])
        
        n_to_remove = safe_split(row['N_to_Remove'])
        if n_to_remove:
            for i in n_to_remove:
                image[image == int(i)] = 0

        mapper = Morph.modules.Mapper()
        cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

        n_label = image
        xenium_annotations[section_id]['n_label'] = n_label

        output_file = os.path.join(section_dir, f"{section_id}_n_label.csv")
        Morph.writers.xenium(output_file, n_label, cells)
        print(f"write {output_file}")

    
        n_mask = (image > 0).astype(int)
        xenium_annotations[section_id]['n_mask'] = n_mask

# 4. Tumor mask
print("\n=== Step 4: Tumor mask ===")

## A. functions
def t_filter(image, mask, element=None):
    image = numpy.multiply(image, mask)
    image = Morph.operators.dilation(image, skimage.morphology.disk(1))
    return image

def t_area_filter(image, lambda_area_closing, lambda_area_opening):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    return skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)

## B. Processing
d = 10
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing Tumor mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    tumor_microregion = safe_str(row['Tumor_microregion'])
    if tumor_microregion is None:
        raise ValueError(f"Tumor_microregion is required for {section_id}")
    G = clean_cluster_set(tumor_microregion)
    if not G:
        raise ValueError(f"Tumor_microregion is empty for {section_id}")
    
    image = Morph.backbone(xenium_annotations[section_id]['data'], 
                           ['xenium', d], 
                           ['total', G], 
                           ['maximum'],  
                           ['custom', t_filter, 1 - xenium_annotations[section_id]['cm_mask'] - xenium_annotations[section_id]['lp_mask'] - xenium_annotations[section_id]['lubnat_mask'], S], 
                           ['binary', tau], 
                           ['custom', t_area_filter, 10000, 30], 
                           ['blob', S])

    t_to_remove = safe_split(row['T_to_remove'])
    if t_to_remove:
        for i in t_to_remove:
            image[image == int(i)] = 0

    mapper = Morph.modules.Mapper()
    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

    t_label = image
    xenium_annotations[section_id]['t_label'] = t_label

    output_file = os.path.join(section_dir, f"{section_id}_t_label.csv")
    Morph.writers.xenium(output_file, t_label, cells)
    print(f"write {output_file}")

    
    t_mask = (image > 0).astype(int)
    xenium_annotations[section_id]['t_mask'] = t_mask
    output_file = os.path.join(section_dir, f"{section_id}_t_mask.csv")
    Morph.writers.xenium(output_file, t_mask, cells)
    print(f"write {output_file}")

# 5. Tumor-Necrosis mask
print("\n=== Step 5: Tumor-Necrosis mask ===")

## A. Functions
def tn_filter(image, t_mask, n_mask, element=None):
    image = numpy.multiply(image, (t_mask + n_mask) > 0)
    image = Morph.operators.dilation(image, skimage.morphology.disk(1))
    return Morph.operators.opening(Morph.operators.closing(image, element), element)

def tn_area_filter(image, lambda_area_closing, lambda_area_opening):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    image = skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)
    return image

def tn_hole_area_filter(image, lambda_area_closing, lambda_area_opening, hole_path):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    image = skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)
    hole_mask = hole_annotation(hole_path, image.shape, d)
    output_image = ((image + hole_mask) > 0).astype(int)
    return output_image

## C. Post-processing
d = 10
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing Tumor-Necrosis mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    necrosis = safe_str(row['Necrosis'])
    tumor_microregion = safe_str(row['Tumor_microregion'])
    if tumor_microregion is None:
        raise ValueError(f"Tumor_microregion is required for {section_id}")
    
    if necrosis is not None: 
        tn = tumor_microregion + ',' + necrosis
        G = clean_cluster_set(tn)
    else:
        G = clean_cluster_set(tumor_microregion)
    if not G:
        raise ValueError(f"Tumor_microregion+Necrosis combination is empty for {section_id}")

    hole_coord = safe_str(row['hole_coordinate'])
    if hole_coord is not None:
        hole_path = hole_coord

        image = Morph.backbone(xenium_annotations[section_id]['data'], 
                           ['xenium', d], 
                           ['total', G], 
                           ['maximum'],  
                           ['custom', tn_filter, xenium_annotations[section_id]['t_mask'], xenium_annotations[section_id]['n_mask'], S], 
                           ['binary', tau], 
                           ['custom', tn_hole_area_filter, 10000, 50, hole_path], 
                           ['blob', S])
    else:
        image = Morph.backbone(xenium_annotations[section_id]['data'], 
                               ['xenium', d], 
                               ['total', G], 
                               ['maximum'],  
                               ['custom', tn_filter, xenium_annotations[section_id]['t_mask'], xenium_annotations[section_id]['n_mask'], S], 
                               ['binary', tau], 
                               ['custom', tn_area_filter, 10000, 50], 
                               ['blob', S])

    tn_to_remove = safe_split(row['TN_to_remove'])
    if tn_to_remove:
        print(tn_to_remove)
        for i in tn_to_remove:
            image[image == int(i)] = 0    

    mapper = Morph.modules.Mapper()
    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

    tn_label = image
    xenium_annotations[section_id]['tn_label'] = tn_label

    output_file = os.path.join(section_dir, f"{section_id}_tn_label.csv")
    Morph.writers.xenium(output_file, tn_label, cells)
    print(f"write {output_file}")

    
    tn_mask = (image > 0).astype(int)
    xenium_annotations[section_id]['tn_mask'] = tn_mask
    output_file = os.path.join(section_dir, f"{section_id}_tn_mask.csv")
    Morph.writers.xenium(output_file, tn_mask, cells)
    print(f"write {output_file}")

# 6. S mask
print("\n=== Step 6: S mask ===")

## A. Functions
def s_filter(image, mask, dilation_n, element=None):
    image = numpy.multiply(image, mask)
    image = Morph.operators.dilation(image, skimage.morphology.disk(dilation_n))
    return Morph.operators.opening(Morph.operators.closing(image, element), element)

def s_area_filter(image, mask, lambda_area_closing, lambda_area_opening):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    
    image = ((image - mask) > 0)
    
    output_image = skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)
    return output_image

## C. Post-Processing
d = 10
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing S (Stroma) mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    gross_tumor_stroma = safe_str(row['Gross_tumor_stroma'])
    if gross_tumor_stroma is None:
        raise ValueError(f"Gross_tumor_stroma is required for {section_id}")
    G = clean_cluster_set(gross_tumor_stroma)
    if not G:
        raise ValueError(f"Gross_tumor_stroma is empty for {section_id}")

    filtered_mask = ((xenium_annotations[section_id]['cm_mask'] + xenium_annotations[section_id]['cl_mask'] + xenium_annotations[section_id]['c_mus_mask'] + xenium_annotations[section_id]['tn_mask'] + xenium_annotations[section_id]['lp_mask'] + xenium_annotations[section_id]['lubnat_mask']) > 0).astype(int)

    dilation_n = 2

    if section_id in ('HT472C1-M1FP1U1'):
        dilation_n = 0 

    image = Morph.backbone(xenium_annotations[section_id]['data'], 
                           ['xenium', d], 
                           ['total', G], 
                           ['maximum'],  
                           ['custom', s_filter, 1 - filtered_mask, dilation_n, S], 
                           ['binary', tau], 
                           ['custom', s_area_filter, filtered_mask, 10000, 300], 
                           ['blob', S])

    s_to_remove = safe_split(row['S_to_remove'])
    if s_to_remove:
        print(s_to_remove)
        for i in s_to_remove:
            image[image == int(i)] = 0 

    mapper = Morph.modules.Mapper()
    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

    s_label = image
    xenium_annotations[section_id]['s_label'] = s_label

    output_file = os.path.join(section_dir, f"{section_id}_s_label.csv")
    Morph.writers.xenium(output_file, s_label, cells)
    print(f"write {output_file}")


    s_mask = (image > 0).astype(int)
    xenium_annotations[section_id]['s_mask'] = s_mask
    output_file = os.path.join(section_dir, f"{section_id}_s_mask.csv")
    Morph.writers.xenium(output_file, s_mask, cells)
    print(f"write {output_file}")

# 7. TNS mask
print("\n=== Step 7: TNS mask ===")

def tns_filter(image, tn_mask, s_mask, element=None):
    image = numpy.multiply(image, (tn_mask + s_mask) > 0)
    image = Morph.operators.dilation(image, skimage.morphology.disk(1))
    return Morph.operators.opening(Morph.operators.closing(image, element), element)

def tns_area_filter(image, lambda_area_closing, lambda_area_opening):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    image = skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)
    return image

def tns_hole_area_filter(image, lambda_area_closing, lambda_area_opening, hole_path):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    image = skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)
    hole_mask = hole_annotation(hole_path, image.shape, d)
    output_image = ((image + hole_mask) > 0).astype(int)
    return output_image

def tns_hole_area_filter2(image, lambda_area_closing, lambda_area_opening, hole_path, remove_cor):
    image = skimage.morphology.remove_small_holes(image, lambda_area_closing)
    image = skimage.morphology.remove_small_objects(image, lambda_area_opening).astype(int)
    hole_mask = hole_annotation(hole_path, image.shape, d)
    remove_mask = hole_annotation(remove_cor, image.shape, d)
    output_image = ((image + hole_mask - remove_mask) > 0).astype(int)
    return output_image

d = 10
S = numpy.ones((3, 3))
tau = 1

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing TNS (Tumor-Necrosis-Stroma) mask for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    necrosis = safe_str(row['Necrosis'])
    tumor_microregion = safe_str(row['Tumor_microregion'])
    gross_tumor_stroma = safe_str(row['Gross_tumor_stroma'])
    if tumor_microregion is None:
        raise ValueError(f"Tumor_microregion is required for {section_id}")
    if gross_tumor_stroma is None:
        raise ValueError(f"Gross_tumor_stroma is required for {section_id}")
    
    if necrosis is not None: 
        tns = tumor_microregion + ',' + necrosis + ',' + gross_tumor_stroma
        G = clean_cluster_set(tns)
    else:
        tns = tumor_microregion + ',' + gross_tumor_stroma
        G = clean_cluster_set(tns)
    if not G:
        raise ValueError(f"Tumor_microregion+Necrosis+Gross_tumor_stroma combination is empty for {section_id}")

    tn_mask = xenium_annotations[section_id]['tn_mask']
    s_mask = xenium_annotations[section_id]['s_mask']
    
    
    tns_hole_coord = safe_str(row['TNS_hole_coordinate'])
    if tns_hole_coord is not None:
        tns_remove_coord = safe_str(row['TNS_remove_coordinate'])
        if tns_remove_coord is not None:
            hole_path = tns_hole_coord
            remove_cor_path = tns_remove_coord
            image = Morph.backbone(xenium_annotations[section_id]['data'], 
                   ['xenium', d], 
                   ['total', G], 
                   ['maximum'],  
                   ['custom', tns_filter, tn_mask, s_mask, S], 
                   ['binary', tau], 
                   ['custom', tns_hole_area_filter2, 1000, 500, hole_path, remove_cor_path], 
                   ['blob', S])
        else:    
            hole_path = tns_hole_coord

            image = Morph.backbone(xenium_annotations[section_id]['data'], 
                               ['xenium', d], 
                               ['total', G], 
                               ['maximum'],  
                               ['custom', tns_filter, tn_mask, s_mask, S], 
                               ['binary', tau], 
                               ['custom', tns_hole_area_filter, 1000, 500, hole_path], 
                               ['blob', S])
        

    else:
        image = Morph.backbone(xenium_annotations[section_id]['data'], 
                               ['xenium', d], 
                               ['total', G], 
                               ['maximum'],  
                               ['custom', tns_filter, tn_mask, s_mask, S], 
                               ['binary', tau], 
                               ['custom', tns_area_filter, 1000, 500], 
                               ['blob', S])

    tns_to_remove = safe_split(row['TNS_to_remove'])
    if tns_to_remove:
        print(tns_to_remove)
        for i in tns_to_remove:
            image[image == int(i)] = 0 

    mapper = Morph.modules.Mapper()
    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

    tns_label = image
    xenium_annotations[section_id]['tns_label'] = tns_label

    output_file = os.path.join(section_dir, f"{section_id}_tns_label.csv")
    Morph.writers.xenium(output_file, tns_label, cells)
    print(f"write {output_file}")

    
    tns_mask = (image > 0).astype(int)
    xenium_annotations[section_id]['tns_mask'] = tns_mask
    output_file = os.path.join(section_dir, f"{section_id}_tns_mask.csv")
    Morph.writers.xenium(output_file, tns_mask, cells)
    print(f"write {output_file}")

# 8. TN Distance
print("\n=== Step 8: TN Distance ===")

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing TN Distance for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    file_path = safe_str(row['Tissue_margin'])
    if file_path is None:
        raise ValueError(f"Tissue_margin is required for {section_id}")
    margin_mask = annotation(file_path, xenium_annotations[section_id]['tn_mask'].shape, d)

    inward_distances = Distance().maximum(1 - (xenium_annotations[section_id]['tn_mask'] > 0).astype(int), tissue=margin_mask).astype(int)
    outward_distances = Distance().maximum(xenium_annotations[section_id]['tn_mask'].astype(int), tissue=margin_mask).astype(int)

    mapper = Morph.modules.Mapper()
    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)    

    inward_file = os.path.join(section_dir, f"{section_id}_tn_inward.csv")
    inward_distances[inward_distances < 0] = 0
    Morph.writers.xenium(inward_file, inward_distances, cells)   

    outward_file = os.path.join(section_dir, f"{section_id}_tn_outward.csv")
    outward_distances[outward_distances < 0] = 0
    Morph.writers.xenium(outward_file, outward_distances, cells)   
    
    print(f"write {section_id}")

# 9. TNS Distance
print("\n=== Step 9: TNS Distance ===")

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing TNS Distance for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    file_path = safe_str(row['Tissue_margin'])
    if file_path is None:
        raise ValueError(f"Tissue_margin is required for {section_id}")
    margin_mask = annotation(file_path, xenium_annotations[section_id]['tns_mask'].shape, d)

    inward_distances = Distance().maximum(1 - (xenium_annotations[section_id]['tns_mask'] > 0).astype(int), tissue=margin_mask).astype(int)
    outward_distances = Distance().maximum(xenium_annotations[section_id]['tns_mask'].astype(int), tissue=margin_mask).astype(int)

    mapper = Morph.modules.Mapper()
    cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)    

    inward_file = os.path.join(section_dir, f"{section_id}_tns_inward.csv")
    inward_distances[inward_distances < 0] = 0
    Morph.writers.xenium(inward_file, inward_distances, cells)   

    outward_file = os.path.join(section_dir, f"{section_id}_tns_outward.csv")
    outward_distances[outward_distances < 0] = 0
    Morph.writers.xenium(outward_file, outward_distances, cells)   
    
    print(f"write {section_id}")

# 10-1. TN outward per TMR
print("\n=== Step 10: TN outward per TMR ===")

d = 10

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing TN outward per TMR for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    file_path = safe_str(row['Tissue_margin'])
    if file_path is None:
        raise ValueError(f"Tissue_margin is required for {section_id}")
    tn_mask_shape = xenium_annotations[section_id]['tn_mask'].shape
    margin_mask = annotation(file_path, tn_mask_shape, d)

    tn_labels = numpy.unique(xenium_annotations[section_id]['tn_label'])
    tn_labels = tn_labels[tn_labels > 0]  # Ignore background (assuming 0 means background)

    for label in tn_labels:
        label_mask = (xenium_annotations[section_id]['tn_label'] == label).astype(int)
        distances_out = Distance().maximum(label_mask).astype(int)
        distances_out = numpy.multiply(distances_out, margin_mask)

        distances_out[distances_out > 11] = 0
        distances_out[distances_out < 0] = 0

        mapper = Morph.modules.Mapper()
        cells = mapper.xenium(xenium_annotations[section_id]['cell'], d)

        os.makedirs(os.path.join(section_dir, 'tmr_outward'), exist_ok=True)
        outward_file = os.path.join(section_dir, 'tmr_outward', f"{section_id}_tn{label}_outward.csv")
        Morph.writers.xenium(outward_file, distances_out, cells)

        print(f"Written: {section_id} - tn{label}_outward.csv")

# 11. Size and Shape
print("\n=== Step 11: Size and Shape ===")

summary_data = []

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    print(f"\n  Processing Size and Shape for {section_id}...")
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    size = Morph.features.Size()
    
    # Extract scalar values from dictionary output
    t_mask_size = list(size.count(xenium_annotations[section_id]['t_mask']).values())[0]
    tn_mask_size = list(size.count(xenium_annotations[section_id]['tn_mask']).values())[0]
    tns_mask_size = list(size.count(xenium_annotations[section_id]['tns_mask']).values())[0]

    # Collect data
    summary_data.append({
        "section_id": section_id,
        "t_mask_size": t_mask_size,
        "tn_mask_size": tn_mask_size,
        "tns_mask_size": tns_mask_size
    })

    print(f"Processed {section_id}")

# Create and save the summary DataFrame
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(os.path.join(OUTDIR, "mask_size_summary.csv"), index=False)

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    size = Morph.features.Size()
    t_size = size.count(xenium_annotations[section_id]['t_label']) 

    t_size_file = os.path.join(section_dir, f"{section_id}_t_size.csv")
    Morph.writers.xenium_dict(t_size_file, t_size)
    
    print(f"write {section_id}")

for index, row in mask_df.iterrows():
    section_id = row['SectionID']
    section_dir = os.path.join(OUTDIR, section_id)
    os.makedirs(section_dir, exist_ok=True)

    size = Morph.features.Size()
    tn_size = size.count(xenium_annotations[section_id]['tn_label']) 
    tns_size = size.count(xenium_annotations[section_id]['tns_label']) 

    tn_size_file = os.path.join(section_dir, f"{section_id}_tn_size.csv")
    Morph.writers.xenium_dict(tn_size_file, tn_size)

    tns_size_file = os.path.join(section_dir, f"{section_id}_tns_size.csv")
    Morph.writers.xenium_dict(tns_size_file, tns_size)
  
    
    print(f"write {section_id}")

# for index, row in mask_df.iterrows():
#     start_time = time.time()

#     section_id = row['SectionID']
#     section_dir = os.path.join(OUTDIR, section_id)
#     os.makedirs(section_dir, exist_ok=True)

#     shape = Morph.features.Shape()
#     tn_shape = shape.roundness(xenium_annotations[section_id]['tn_label']) 

#     tn_shape_file = os.path.join(section_dir, f"{section_id}_tn_shape.csv")
#     Morph.writers.xenium_dict(tn_shape_file, tn_shape)

#     elapsed = time.time() - start_time
#     print(f"Finished processing {section_id} in {elapsed:.2f} seconds")

# 12. save dictionary
print("\n=== Step 12: Saving annotations ===")

morph_annotation_path = os.path.join(OUTDIR, "xenium_morph_annotations.pkl")

with open(morph_annotation_path, "wb") as f:
    pickle.dump(xenium_annotations, f)

print(f"Saved annotations to {morph_annotation_path}")
print("\n=== Processing Complete ===")

