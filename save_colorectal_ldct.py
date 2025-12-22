from algorithm.sim_functions import *
from algorithm.utils_functions import *
import os
import numpy as np
import pydicom
from pydicom.uid import generate_uid
import matplotlib.pyplot as plt

def main():
    
    base_col_path = '/home/data/Datasets/Colorectal-Liver-Metastases/'
    patients_with_mas, list_patients_with_mas = col_mas_av(base_col_path)

    # Carpetas destino para guardar las simulaciones
    save_dir_ld2_p = "/home/data/Datasets/Colorectal-Liver-Metastases//crlm_ld_2_poisson"
    save_dir_ld4_p = "/home/data/Datasets/Colorectal-Liver-Metastases//crlm_ld_4_poisson"
    save_dir_ld8_p = "/home/data/Datasets/Colorectal-Liver-Metastases//crlm_ld_8_poisson"
    save_dir_ld2_pg = "/home/data/Datasets/Colorectal-Liver-Metastases//crlm_ld_2_gauss_poisson"
    save_dir_ld4_pg = "/home/data/Datasets/Colorectal-Liver-Metastases//crlm_ld_4_gauss_poisson"
    save_dir_ld8_pg = "/home/data/Datasets/Colorectal-Liver-Metastases//crlm_ld_8_gauss_poisson"

    os.makedirs(save_dir_ld2_p, exist_ok=True)
    os.makedirs(save_dir_ld4_p, exist_ok=True)
    os.makedirs(save_dir_ld8_p, exist_ok=True)
    os.makedirs(save_dir_ld2_pg, exist_ok=True)
    os.makedirs(save_dir_ld4_pg, exist_ok=True)
    os.makedirs(save_dir_ld8_pg, exist_ok=True)

    # ----------------------------------------------------------
    # Función auxiliar para guardar un DICOM con píxeles nuevos
    # ----------------------------------------------------------
    def save_lowdose_dicom(original_path, low_img, out_path, new_mas):
        ds = pydicom.dcmread(original_path)

        # Normalizar low_img al rango original del DICOM
        low_norm = (low_img - np.min(low_img)) / (np.max(low_img) - np.min(low_img))
        low_norm = (low_norm * 4095).astype(np.uint16)

        ds.PixelData = low_norm.tobytes()

        # -------------------------------------------------
        # Actualizar mAs si existe el tag Exposure (0018,1152)
        # -------------------------------------------------
        if hasattr(ds, "Exposure"):
            ds.Exposure = float(new_mas)
        else:
            # Crear el tag manualmente si no existe
            from pydicom.tag import Tag
            ds.add_new(Tag(0x0018, 0x1152), 'DS', float(new_mas))

        # Si quieres log: print("Nuevo mAs guardado:", new_mas)

        # Mantener UIDs o generar nuevos
        ds.SOPInstanceUID = generate_uid()

        # Guardar
        ds.save_as(out_path)



    # ==========================================================
    #     PIPELINE PRINCIPAL
    # ==========================================================

    list_h_p = []
    list_l_p = []
    list_h_pg = []
    list_l_pg = []

    for pat_name, pat_mas in patients_with_mas:

        pat_path = os.path.join(base_col_path, pat_name)
        list_paths = sorted(os.listdir(pat_path))
        dcm_path = os.path.join(pat_path, list_paths[0])

        print(f"\n=== Procesando paciente {pat_name} (mAs original: {pat_mas}) ===")

        list_seg_dcm = os.listdir(dcm_path)
        for p in list_seg_dcm:
            if 'Segmentation' in p:
                continue

            high_paths = os.path.join(dcm_path, p)
            list_dcms = sorted(os.listdir(high_paths))

            # Vista previa del primer slice
            #dcm_0 = os.path.join(high_paths, list_dcms[0])
            #array_dcm = pydicom.dcmread(dcm_0).pixel_array
            #plt.title(f'Patient: {pat_name} - MAS: {pat_mas}')
            #plt.imshow(array_dcm, cmap='gray')
            #plt.show()

            # Subcarpetas por paciente en las simulaciones
            pat_ld2_path_p = os.path.join(save_dir_ld2_p, pat_name)
            pat_ld4_path_p = os.path.join(save_dir_ld4_p, pat_name)
            pat_ld8_path_p = os.path.join(save_dir_ld8_p, pat_name)
            pat_ld2_path_pg = os.path.join(save_dir_ld2_pg, pat_name)
            pat_ld4_path_pg = os.path.join(save_dir_ld4_pg, pat_name)
            pat_ld8_path_pg = os.path.join(save_dir_ld8_pg, pat_name)   
            os.makedirs(pat_ld2_path_p, exist_ok=True)
            os.makedirs(pat_ld4_path_p, exist_ok=True)
            os.makedirs(pat_ld8_path_p, exist_ok=True)
            os.makedirs(pat_ld2_path_pg, exist_ok=True)
            os.makedirs(pat_ld4_path_pg, exist_ok=True)
            os.makedirs(pat_ld8_path_pg, exist_ok=True)

            # Listas temporales
            list_h_dcms_p1, list_l_dcms_p1 = [], []
            list_h_dcms_p2, list_l_dcms_p2 = [], []
            list_h_dcms_p3, list_l_dcms_p3 = [], []
            list_h_dcms_pg1, list_l_dcms_pg1 = [], []
            list_h_dcms_pg2, list_l_dcms_pg2 = [], []
            list_h_dcms_pg3, list_l_dcms_pg3 = [], []

            for dc in list_dcms:

                high_path = os.path.join(high_paths, dc)

                mas_des1 = pat_mas // 2
                mas_des2 = pat_mas // 4
                mas_des3 = pat_mas // 8

                # ---- Simulaciones ----
                h_p1, low_p1, _ = apply_ld_sinogram_full(high_path, mas_des1, print_logs=False, mod=False, gaussian_noise=False)
                h_p2, low_p2, _ = apply_ld_sinogram_full(high_path, mas_des2, print_logs=False, mod=False, gaussian_noise=False)
                h_p3, low_p3, _ = apply_ld_sinogram_full(high_path, mas_des3, print_logs=False, mod=False, gaussian_noise=False)
                h_pg1, low_pg1, _ = apply_ld_sinogram_full(high_path, mas_des1, print_logs=False, mod=False, gaussian_noise=True)
                h_pg2, low_pg2, _ = apply_ld_sinogram_full(high_path, mas_des2, print_logs=False, mod=False, gaussian_noise=True)
                h_pg3, low_pg3, _ = apply_ld_sinogram_full(high_path, mas_des3, print_logs=False, mod=False, gaussian_noise=True)
                # ---- Guardado DICOM ----
                out_p1 = os.path.join(pat_ld2_path_p, f"LD2_{dc}")
                out_p2 = os.path.join(pat_ld4_path_p, f"LD4_{dc}")
                out_p3 = os.path.join(pat_ld8_path_p, f"LD8_{dc}")
                out_pg1 = os.path.join(pat_ld2_path_pg, f"LD2_{dc}")
                out_pg2 = os.path.join(pat_ld4_path_pg, f"LD4_{dc}")
                out_pg3 = os.path.join(pat_ld8_path_pg, f"LD8_{dc}")
                save_lowdose_dicom(high_path, low_p1, out_p1, mas_des1)
                save_lowdose_dicom(high_path, low_p2, out_p2, mas_des2)
                save_lowdose_dicom(high_path, low_p3, out_p3, mas_des3)
                save_lowdose_dicom(high_path, low_pg1, out_pg1, mas_des1)
                save_lowdose_dicom(high_path, low_pg2, out_pg2, mas_des2)
                save_lowdose_dicom(high_path, low_pg3, out_pg3, mas_des3)


                # Guardar temporal en memoria
                list_h_dcms_p1.append(h_p1)
                list_l_dcms_p1.append(low_p1)
                list_h_dcms_p2.append(h_p2)
                list_l_dcms_p2.append(low_p2)
                list_h_dcms_p3.append(h_p3)
                list_l_dcms_p3.append(low_p3)
                list_h_dcms_pg1.append(h_pg1)
                list_l_dcms_pg1.append(low_pg1)
                list_h_dcms_pg2.append(h_pg2)
                list_l_dcms_pg2.append(low_pg2)
                list_h_dcms_pg3.append(h_pg3)
                list_l_dcms_pg3.append(low_pg3)

            # ==========================================================
            #     VISUALIZAR UN SLICE CENTRAL (alta / baja dosis)
            # ==========================================================
            idx_mid = len(list_h_dcms_p1) // 2

            plt.figure(figsize=(12, 6))
            plt.suptitle(f"Paciente: {pat_name} – Comparación de dosis (Gauss)", fontsize=14)

            plt.subplot(1, 4, 1)
            plt.imshow(list_h_dcms_p1[idx_mid], cmap='gray')
            plt.title(f"Alta dosis\n({pat_mas} mAs)")
            plt.axis('off')

            plt.subplot(1, 4, 2)
            plt.imshow(list_l_dcms_p1[idx_mid], cmap='gray')
            plt.title(f"½ mAs\n({mas_des1} mAs)")
            plt.axis('off')

            plt.subplot(1, 4, 3)
            plt.imshow(list_l_dcms_p2[idx_mid], cmap='gray')
            plt.title(f"¼ mAs\n({mas_des2} mAs)")
            plt.axis('off')

            plt.subplot(1, 4, 4)
            plt.imshow(list_l_dcms_p3[idx_mid], cmap='gray')
            plt.title(f"⅛ mAs\n({mas_des3} mAs)")
            plt.axis('off')

            plt.show()

            idx_mid = len(list_h_dcms_pg1) // 2

            plt.figure(figsize=(12, 6))
            plt.suptitle(f"Paciente: {pat_name} – Comparación de dosis (Gauss+Poisson)", fontsize=14)

            plt.subplot(1, 4, 1)
            plt.imshow(list_h_dcms_pg1[idx_mid], cmap='gray')
            plt.title(f"Alta dosis\n({pat_mas} mAs)")
            plt.axis('off')

            plt.subplot(1, 4, 2)
            plt.imshow(list_l_dcms_pg1[idx_mid], cmap='gray')
            plt.title(f"½ mAs\n({mas_des1} mAs)")
            plt.axis('off')

            plt.subplot(1, 4, 3)
            plt.imshow(list_l_dcms_pg2[idx_mid], cmap='gray')
            plt.title(f"¼ mAs\n({mas_des2} mAs)")
            plt.axis('off')

            plt.subplot(1, 4, 4)
            plt.imshow(list_l_dcms_pg3[idx_mid], cmap='gray')
            plt.title(f"⅛ mAs\n({mas_des3} mAs)")
            plt.axis('off')

            plt.show()

            # Guardar para análisis
            list_h_p.append(np.array(list_h_dcms_p1))
            list_l_p.append(np.array(list_l_dcms_p1))
            list_h_pg.append(np.array(list_h_dcms_pg1))
            list_l_pg.append(np.array(list_l_dcms_pg1))

    print("\n=== PROCESO COMPLETADO ===")

if __name__ == "__main__":
    main()
