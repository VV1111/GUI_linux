# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
# """
# Simple Streamlit-based Medical Image Viewer

# - Scan a data root for NIfTI (.nii / .nii.gz) files
# - Group files by case
# - Select case & slice in sidebar
# - Display 2D axial slice
# - (Optional) Overlay segmentation mask if available

# 适合作为 TopoVASC / CTA/MRA 查看工具的原型。
# """

# import os
# from typing import Dict, List, Optional, Tuple

# import numpy as np
# import nibabel as nib
# import streamlit as st
# import matplotlib.pyplot as plt

# # ==============================
# # 1. 配置：数据根目录
# # ==============================
# # TODO: 修改为你自己的数据路径，例如：
# # /home/zhiwei/research/segmentation_all/TopoVASC/dataset/nii
# DATA_ROOT = "/home/zhiwei/research/data/ESUS"


# # ==============================
# # 2. 工具函数：扫描病例
# # ==============================
# def scan_cases(data_root: str) -> Dict[str, Dict[str, str]]:
#     """
#     扫描数据目录，返回:
#     {
#       case_id: {
#           "image": "/path/to/image.nii.gz",
#           "label": "/path/to/label.nii.gz" (optional)
#       },
#       ...
#     }

#     这里给出一个简单的规则示例：
#     - 假设图像文件名类似: case001.nii.gz 或 case001_img.nii.gz
#     - 标签文件名类似: case001_seg.nii.gz 或 case001_label.nii.gz
#     你可以根据自己的命名规范进行修改。
#     """
#     cases: Dict[str, Dict[str, str]] = {}

#     if not os.path.isdir(data_root):
#         return cases

#     for fname in os.listdir(data_root):
#         if not (fname.endswith(".nii") or fname.endswith(".nii.gz")):
#             continue

#         fpath = os.path.join(data_root, fname)
#         basename = fname.replace(".nii.gz", "").replace(".nii", "")

#         # 简单约定：带 "seg" 或 "label" 视为标签
#         is_label = any(k in basename.lower() for k in ["seg", "label", "mask"])

#         # 去掉后缀中的 _seg / _label 等，作为 case_id
#         case_id = basename
#         for key in ["_seg", "_label", "_mask", "-seg", "-label"]:
#             case_id = case_id.replace(key, "")

#         if case_id not in cases:
#             cases[case_id] = {}

#         if is_label:
#             cases[case_id]["label"] = fpath
#         else:
#             cases[case_id]["image"] = fpath

#     # 只保留至少有 image 的 case
#     cases = {k: v for k, v in cases.items() if "image" in v}
#     return cases


# # ==============================
# # 3. 加载 NIfTI
# # ==============================
# @st.cache_data(show_spinner=True)
# def load_nifti(path: str) -> Tuple[np.ndarray, Tuple[float, float, float]]:
#     """
#     读取 NIfTI 文件，返回:
#     - data: numpy array (float32)
#     - spacing: voxel spacing (dx, dy, dz)
#     """
#     img = nib.load(path)
#     data = img.get_fdata().astype(np.float32)
#     header = img.header
#     zooms = header.get_zooms()[:3]
#     return data, zooms


# # ==============================
# # 4. 主界面逻辑
# # ==============================
# def main():
#     st.set_page_config(page_title="Streamlit Medical Image Viewer",
#                        layout="wide")

#     st.title("🧠 Streamlit Medical Image Viewer")
#     st.write(
#         "A simple NIfTI viewer running on a headless Linux server. "
#         "适合作为 TopoVASC / CTA / MRA 的 Web 查看原型。"
#     )

#     # ---- Sidebar: 数据路径 & 病例选择 ----
#     st.sidebar.header("Data Settings / 数据设置")

#     data_root = st.sidebar.text_input("Data root (NIfTI folder)",
#                                       value=DATA_ROOT)
#     st.sidebar.write(f"Current data root: `{data_root}`")

#     cases = scan_cases(data_root)
#     if not cases:
#         st.warning("No NIfTI files found. 请检查 DATA_ROOT 路径或文件命名。")
#         st.stop()

#     case_ids = sorted(list(cases.keys()))
#     case_id = st.sidebar.selectbox("Select case / 选择病例", case_ids)

#     case_info = cases[case_id]
#     img_path = case_info["image"]
#     has_label = "label" in case_info

#     st.sidebar.markdown(f"**Image path:** `{os.path.basename(img_path)}`")
#     if has_label:
#         st.sidebar.markdown(
#             f"**Label path:** `{os.path.basename(case_info['label'])}`"
#         )

#     # 选择显示模式：仅图像 / 图像+标签
#     view_mode = "Image only"
#     if has_label:
#         view_mode = st.sidebar.selectbox(
#             "View mode / 显示模式",
#             ["Image only", "Image + Label overlay"]
#         )

#     # ---- 加载数据 ----
#     img_data, img_spacing = load_nifti(img_path)

#     # 默认按照 axial 方向显示（第三个维度）
#     z_max = img_data.shape[2] - 1
#     z_idx = st.sidebar.slider(
#         "Axial slice index / 轴向切片编号",
#         min_value=0,
#         max_value=z_max,
#         value=z_max // 2,
#         step=1,
#     )

#     # 简单 normalization
#     vmin, vmax = np.percentile(img_data, [1, 99])
#     slice_img = img_data[:, :, z_idx]
#     slice_img = np.clip(slice_img, vmin, vmax)

#     # 读取标签（如果需要 overlay）
#     if has_label and view_mode == "Image + Label overlay":
#         label_data, _ = load_nifti(case_info["label"])
#         if label_data.shape != img_data.shape:
#             st.error(
#                 f"Label shape {label_data.shape} != image shape {img_data.shape}. "
#                 "请检查是否已经对齐或配准。"
#             )
#             label_slice = None
#         else:
#             label_slice = label_data[:, :, z_idx]
#     else:
#         label_slice = None

#     # ---- 主显示区 ----
#     col1, col2 = st.columns([3, 2])

#     with col1:
#         st.subheader(f"Case: {case_id}  |  Axial slice: {z_idx}")
#         fig, ax = plt.subplots(figsize=(5, 5))
#         ax.imshow(slice_img.T, cmap="gray", origin="lower")

#         if label_slice is not None:
#             # 将非零标签叠加着色，这里简单用红色透明覆盖
#             mask = label_slice.T > 0
#             ax.imshow(
#                 np.ma.masked_where(~mask, mask),
#                 cmap="autumn",
#                 alpha=0.4,
#                 origin="lower"
#             )

#         ax.set_axis_off()
#         st.pyplot(fig)

#     with col2:
#         st.subheader("Metadata / 元数据")
#         st.write(f"**Image shape:** {img_data.shape}")
#         st.write(f"**Voxel spacing (dx, dy, dz):** {img_spacing}")

#         if label_slice is not None:
#             unique_labels = np.unique(label_slice)
#             st.write(f"**Unique label values in this slice:** {unique_labels}")

#         st.markdown("---")
#         st.markdown("### Notes / 备注")
#         st.write(
#             "- 可以在此基础上扩展：\n"
#             "  - 三视图（Axial/Coronal/Sagittal）\n"
#             "  - TopoVASC 模型推理按钮\n"
#             "  - Centerline / 拓扑可视化\n"
#             "  - 保存标注结果为 NIfTI / YAML\n"
#         )


# if __name__ == "__main__":
#     main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Streamlit Medical Image Viewer (NIfTI)

This application provides:
- Automatic scanning of NIfTI (.nii / .nii.gz) files in a directory
- Case selection based on file grouping
- Axial slice browsing with a slider
- Optional segmentation overlay (if a label file exists)
- Image metadata display

This serves as a minimal Web-based medical image viewer prototype
that runs on a headless Linux server and is accessible through a browser.
"""

import os
from typing import Dict, Tuple
import numpy as np
import nibabel as nib
import streamlit as st
import matplotlib.pyplot as plt


# ================================================================
# 1. Configuration: Root directory containing NIfTI files
# ================================================================
# TODO: Change this to the directory where your NIfTI files are stored.
DATA_ROOT = "/path/to/your/nifti_root"


# ================================================================
# 2. Scan directory and group files into cases
# ================================================================
def scan_cases(data_root: str) -> Dict[str, Dict[str, str]]:
    """
    Scan the given directory for NIfTI files and group them by case ID.

    Expected naming conventions (flexible and can be customized):
    - Image files:   case001.nii.gz, case001_img.nii.gz, etc.
    - Label files:   case001_seg.nii.gz, case001_label.nii.gz, etc.

    Returns a dictionary of the form:
    {
      "case001": {
          "image": "/path/to/case001.nii.gz",
          "label": "/path/to/case001_seg.nii.gz" (optional)
      },
      ...
    }
    """
    cases: Dict[str, Dict[str, str]] = {}

    if not os.path.isdir(data_root):
        return cases

    for fname in os.listdir(data_root):
        if not (fname.endswith(".nii") or fname.endswith(".nii.gz")):
            continue

        fpath = os.path.join(data_root, fname)
        basename = fname.replace(".nii.gz", "").replace(".nii", "")

        # Heuristic: treat files containing seg/label/mask as labels
        is_label = any(key in basename.lower() for key in ["seg", "label", "mask"])

        # Normalize case ID by removing suffixes such as _seg / _label
        case_id = basename
        for key in ["_seg", "_label", "_mask", "-seg", "-label"]:
            case_id = case_id.replace(key, "")

        if case_id not in cases:
            cases[case_id] = {}

        if is_label:
            cases[case_id]["label"] = fpath
        else:
            cases[case_id]["image"] = fpath

    # Keep only cases that have an image file
    cases = {cid: info for cid, info in cases.items() if "image" in info}
    return cases


# ================================================================
# 3. Load NIfTI volume
# ================================================================
@st.cache_data(show_spinner=True)
def load_nifti(path: str) -> Tuple[np.ndarray, Tuple[float, float, float]]:
    """
    Load a NIfTI file and return:
    - data: numpy array of shape (H, W, D)
    - spacing: voxel spacing (dx, dy, dz)
    """
    img = nib.load(path)
    data = img.get_fdata().astype(np.float32)
    spacing = img.header.get_zooms()[:3]
    return data, spacing


# ================================================================
# 4. Main Streamlit application
# ================================================================
def main():
    st.set_page_config(page_title="Streamlit Medical Image Viewer", layout="wide")

    st.title("🩺 Medical Image Viewer (NIfTI)")
    st.write(
        "A lightweight Web-based viewer for medical images, running on a "
        "headless Linux server and accessible through a browser."
    )

    # ---------- Sidebar: Data root & case selection ----------
    st.sidebar.header("Data Settings")

    data_root = st.sidebar.text_input("NIfTI Directory", value=DATA_ROOT)
    st.sidebar.write(f"Current directory: `{data_root}`")

    # Scan cases
    cases = scan_cases(data_root)
    if not cases:
        st.warning("No NIfTI files found in the specified directory.")
        st.stop()

    case_ids = sorted(list(cases.keys()))
    case_id = st.sidebar.selectbox("Select Case", case_ids)

    case_info = cases[case_id]
    img_path = case_info["image"]
    has_label = "label" in case_info

    st.sidebar.markdown(f"**Image file:** `{os.path.basename(img_path)}`")
    if has_label:
        st.sidebar.markdown(f"**Label file:** `{os.path.basename(case_info['label'])}`")

    # View mode options
    view_mode = "Image only"
    if has_label:
        view_mode = st.sidebar.selectbox(
            "View Mode",
            ["Image only", "Image + Segmentation Overlay"]
        )

    # ---------- Load image (and label) ----------
    img_data, img_spacing = load_nifti(img_path)

    z_max = img_data.shape[2] - 1
    slice_idx = st.sidebar.slider(
        "Axial Slice Index",
        min_value=0,
        max_value=z_max,
        value=z_max // 2,
        step=1,
    )

    # Slice extraction and normalization
    slice_img = img_data[:, :, slice_idx]
    vmin, vmax = np.percentile(slice_img, [1, 99])
    slice_img = np.clip(slice_img, vmin, vmax)

    # Load segmentation slice if needed
    if has_label and view_mode == "Image + Segmentation Overlay":
        label_data, _ = load_nifti(case_info["label"])
        if label_data.shape != img_data.shape:
            st.error(
                f"Label shape {label_data.shape} does not match image shape {img_data.shape}. "
                "Please verify alignment or preprocessing."
            )
            label_slice = None
        else:
            label_slice = label_data[:, :, slice_idx]
    else:
        label_slice = None

    # ---------- Layout: image on left, metadata on right ----------
    col1, col2 = st.columns([3, 2])

    with col1:
        st.subheader(f"Case: {case_id}  |  Slice: {slice_idx}")

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.imshow(slice_img.T, cmap="gray", origin="lower")

        if label_slice is not None:
            mask = label_slice.T > 0
            ax.imshow(
                np.ma.masked_where(~mask, mask),
                cmap="autumn",
                alpha=0.4,
                origin="lower",
            )

        ax.set_axis_off()
        st.pyplot(fig)

    with col2:
        st.subheader("Metadata")
        st.write(f"**Image shape:** {img_data.shape}")
        st.write(f"**Voxel spacing (dx, dy, dz):** {img_spacing}")

        if label_slice is not None:
            unique_vals = np.unique(label_slice)
            st.write(f"**Unique label values (this slice):** {unique_vals}")

        st.markdown("---")
        st.markdown("### Notes")
        st.write(
            "- You may extend this viewer with additional functionality:\n"
            "  - Multiplanar views (axial, coronal, sagittal)\n"
            "  - Window/level controls\n"
            "  - Basic annotation tools\n"
            "  - Overlay of model predictions\n"
            "  - Export of images or annotations\n"
        )


if __name__ == "__main__":
    main()
