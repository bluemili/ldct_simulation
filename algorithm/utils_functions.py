import numpy as np
import matplotlib.pyplot as plt
import SimpleITK as sitk
import sys
import os
import pydicom

def print_values(array, name="Array"):
    print(f"------------------Values of {name}----------------------")
    print(f"Max value of {name}:", np.max(array))
    print(f"Min value of {name}:", np.min(array))
    print(f"Mean value of {name}:", np.mean(array))
    print(f"Std value of {name}:", np.std(array))

def histogram_values(array, name="Array", bins=100, log_scale=False):
    """
    Visualiza un histograma con mejor formato y control sobre los parámetros.

    Parámetros:
    -----------
    array : np.ndarray
        Datos numéricos a graficar.
    name : str
        Nombre descriptivo del conjunto de datos (para el título y etiquetas).
    bins : int o secuencia
        Número de bins o bordes del histograma.
    log_scale : bool
        Si es True, usa escala logarítmica en el eje Y.
    """
    # Asegurar que el array sea plano
    array = np.asarray(array).ravel()
    
    # Crear histograma
    hist, bin_edges = np.histogram(array, bins=bins)

    # Crear figura
    plt.figure(figsize=(8, 5))
    plt.bar(
        bin_edges[:-1],
        hist,
        width=np.diff(bin_edges),
        color='steelblue',
        edgecolor='black',
        alpha=0.7
    )
    
    # Etiquetas y formato
    plt.title(f"Distribución de valores - {name}", fontsize=14)
    plt.xlabel("Valor de píxel", fontsize=12)
    plt.ylabel("Frecuencia", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    
    if log_scale:
        plt.yscale('log')
        plt.ylabel("Frecuencia (escala log)")
    
    plt.tight_layout()
    plt.show()


def affine_from_dicom(path_dcms):
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(path_dcms)
    reader.SetFileNames(dicom_names)
    image = reader.Execute()
    origin = np.array(image.GetOrigin())  # Physical coordinates of the image origin
    spacing = np.array(image.GetSpacing()) # Physical spacing between pixels/voxels
    direction = np.array(image.GetDirection()).reshape(3, 3) # Direction cosines of the image axes

    spacing[0] = spacing[0] * 2
    spacing[1] = spacing[1] * 2
    
    # Create an identity 4x4 matrix
    affine_matrix = np.eye(4)

    # Apply scaling and rotation/orientation
    # The direction matrix already incorporates the spacing implicitly when transforming physical points
    # However, for the affine matrix that transforms voxel coordinates to physical coordinates,
    # you need to explicitly multiply the direction by the spacing.
    affine_matrix[:3, :3] = direction @ np.diag(spacing)

    # Set the translation component (origin)
    affine_matrix[:3, 3] = origin

    return affine_matrix


def calculate_snr_cnr_phan(ev_slice):
    #ROIs (coordenadas [y1:y2, x1:x2])
    roi_tejido = ev_slice[342:355, 300:313] 
    roi_tejido2 = ev_slice[342:355, 189:202]
    mu_tejido, sigma_tejido = np.mean(roi_tejido), np.std(roi_tejido)
    mu_tejido2, sigma_tejido2 = np.mean(roi_tejido2), np.std(roi_tejido2)
    snr_tejido1 = mu_tejido / sigma_tejido
    snr_tejido2 = mu_tejido2 / sigma_tejido2
    cnr = (mu_tejido - mu_tejido2) / ((sigma_tejido + sigma_tejido2) / 2)
    return snr_tejido1, snr_tejido2, cnr, sigma_tejido


def get_pitch(ds):
    """
    Obtiene pitch desde la metadata.
    """
    tags = [
        (0x0018, 0x9311),  # SpiralPitchFactor (más confiable)
        (0x0018, 0x9303),  # TableSpeed
        (0x0018, 0x9304),  # TableFeedPerRotation
    ]
    
    # Spiral pitch factor (Siemens)
    if (0x0018, 0x9311) in ds:
        try:
            return float(ds[(0x0018, 0x9311)].value)
        except:
            pass

    # Intentar reconstruir desde table speed / feed
    try:
        if (0x0018, 0x9303) in ds and (0x0018, 0x9304) in ds:
            speed = float(ds[(0x0018, 0x9303)].value)
            feed = float(ds[(0x0018, 0x9304)].value)
            if feed > 0:
                return speed / feed
    except:
        pass
    
    return None

def get_mAs_from_dicom(ds):
    """
    Intenta obtener el mAs desde diferentes tags estándar y específicos de CT.
    Devuelve el valor encontrado o None si no existe.
    """

    pitch = get_pitch(ds)
    if hasattr(ds, "XRayTubeCurrent"):
        # Si existe y además tenemos el tiempo de exposición -> calcular mAs
        if hasattr(ds, "ExposureTime"):
            # ExposureTime normalmente en ms → convertir a segundos
            mA = float(ds.XRayTubeCurrent)
            t_ms = float(ds.ExposureTime)
            mAs = (mA * (t_ms / 1000.0)) / pitch if pitch else (mA * (t_ms / 1000.0))
            return mAs
    
    possible_tags = [
        (0x0018, 0x1153), # Exposure in mAs
        (0x0018, 0x9332), # Exposure Modulation Strength (not really mAs, avoid)
        (0x0040, 0x030E), # Exposure Modifiers (string, ignore)
    ]
    
    for tag in possible_tags:
        if tag in ds:
            try:
                val = float(ds[tag].value)
                return val
            except:
                pass
    
    # Si no se encontró nada confiable
    return None


def average_mAs_from_dicoms(folder):
    """
    Recorre todos los archivos DICOM de la carpeta y calcula el promedio de mAs.
    Si no encuentra ningún valor, devuelve None.
    """
    mAs_values = []

    for file in os.listdir(folder):
        path = os.path.join(folder, file)
        # Filtrar archivos no dicom (aproximado)
        if not os.path.isfile(path):
            continue

        try:
            ds = pydicom.dcmread(path, stop_before_pixels=True)
        except:
            continue  # saltar archivos corruptos

        mAs = get_mAs_from_dicom(ds)
        if mAs is not None and mAs > 0:
            mAs_values.append(mAs)

    if len(mAs_values) == 0:
        return None
    
    return sum(mAs_values) / len(mAs_values)



def col_mas_av(colorectal_path):
    list_patients = sorted(os.listdir(colorectal_path))

    patients_with_mas = []
    list_patients_with_mas = []
    for patient in list_patients:
        #check if is a directory
        if os.path.isdir(os.path.join(colorectal_path, patient)) == False:
            continue
        patient_path = os.path.join(colorectal_path, patient)
        list_dirs = sorted(os.listdir(patient_path))
        list_paths_seg_na = sorted(os.listdir(os.path.join(patient_path, list_dirs[0])))

        #select dir with NA in the name
        for dire in list_paths_seg_na:
            #agregar que sea una carpeta y no un archivo
            if not "Segmentation" in dire:
                slices_path = os.path.join(patient_path, list_dirs[0], dire)
                mas = average_mAs_from_dicoms(slices_path)
                print(f"Patient: {patient}, mAs average: {mas}")
                if mas is None:
                    print("No mAs value found.")
                    pass
                else:
                    patients_with_mas.append((patient, mas))
                    list_patients_with_mas.append(patient)

    return patients_with_mas, list_patients_with_mas