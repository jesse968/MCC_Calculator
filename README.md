MCC_Calculator - OCT Image Processing Toolkit for Mouse Respiratory Analysis
A comprehensive MATLAB-based software suite for quantitative analysis of Optical Coherence Tomography (OCT) images of mouse respiratory epithelium. This toolkit enables researchers to measure critical biological indicators of mucociliary function, including ciliary beat frequency, mucus transport velocity, and airway surface liquid thickness.

Project Overview
This repository contains three specialized tools designed for sequential analysis of mouse respiratory OCT images:
1. MCC_Calculator_4 - Core preprocessing, segmentation, and functional analysis
2. MCC_Calculator_REG - Image sequence filtering and elastic registration
3. MCC_Calculator_ASL - Airway Surface Liquid (ASL) thickness quantification

These tools work together to transform raw OCT bitmap sequences into quantitative metrics for respiratory research.

Tool 1: MCC_Calculator_4 - Core Analysis Platform
Primary Functions
1. Image Sequence Management
- Multi-format BMP image loading and visualization
- Interactive playback controls with frame navigation
- Multiple processing stage management (Original, Cropped, Filtered, Registered)
- Batch export capabilities

2. Preprocessing Pipeline
- Frame Range Selection: Define temporal analysis windows
- ROI Cropping: Extract regions using rectangular or polygonal selections
- Automatic Frame Filtering: Remove low-quality frames via correlation-based thresholding
- Image Registration: Rigid alignment (affine/similarity/projective transformations)
- Auto-Reference Detection: Identify optimal reference frames using correlation analysis

3. Ciliary Beat Frequency (CBF) Analysis
- Dual-region comparison (ciliary vs. background regions)
- Frequency domain analysis with optional smoothing
- Background motion artifact correction
- Visual frequency spectra plotting with peak identification
- Results logging and export

4. Mucociliary Transport (MCT) Analysis
- DFT-based motion tracking between consecutive frames
- Horizontal and vertical velocity component calculation
- Flow direction analysis
- Unit conversion (pixels → μm/min)
- Cumulative displacement visualization

5. Image Segmentation & Flattening
- Interactive tissue boundary delineation with polylines
- Automated layer segmentation with adjustable parameters
- Curved surface flattening to planar representations
- Manual segmentation adjustment
= Layer coordinate export for further analysis

6. Measurement Calibration
- Spatial calibration (pixels to micrometers)
- Temporal settings (frames per second)
- Unit toggling (pixels ↔ micrometers)

Technical Implementation
- Algorithms: Image correlation, DFT registration, signal processing, image warping
- User Interface: Tab-based organization, interactive ROIs, multi-axes display
- Data Management: State preservation, structured result packaging, logging integration

Tool 2: MCC_Calculator_REG - Quality Control & Registration
Core Functionality
1. Multi-Set Batch Processing
- Hierarchical tree structure for managing multiple datasets
- Context menu operations (delete, expand, collapse)
- Batch operations across multiple experimental conditions
- Visual processing state indicators

2. Advanced Frame Filtering
- Correlation-based quality assessment using normalized cross-correlation
- Adaptive thresholding (user-defined or auto-calculated)
- Visual correlation score plotting
- Frame reduction statistics display

3. Elastic Image Registration
- Multi-scale Demons registration algorithm (imregdemons)
- Pyramid processing for improved accuracy
- Parallel computing for batch acceleration
- Adjustable pyramid levels and smoothing parameters
- Tissue movement and breathing artifact compensation

4. Intelligent Data Organization
- Automatic detection of BMP images in folders/subfolders
- Natural sorting of image sequences
- Separate state preservation for original/filtered/registered images
- Easy switching between processing stages

Technical Features
- Parallel Processing: MATLAB Parallel Computing Toolbox integration
- Quality Control: Correlation plots and frame statistics
- Professional Interface: Dark theme with icon-based navigation
- Scalable Architecture: Efficient handling of large datasets

Tool 3: MCC_Calculator_ASL - Airway Surface Liquid Analysis
Primary Capabilities
1. Multi-Layer Segmentation
- Interactive drawing of surface and basement membrane boundaries
- Automated layer extraction using interpolation algorithms
- Real-time segmentation visualization
- "Save & Next" workflow for batch processing

2. Image Enhancement Pipeline
- Image Flattening: Alignment to reference plane (typically 400-pixel line)
- Adaptive Histogram Equalization: Local contrast enhancement
- Frame Averaging: Temporal smoothing across sequences
- Standard Deviation Mapping: Variability visualization

3. ASL Thickness Measurement
- Epithelial layer top surface detection using gradient-based algorithms
- ASL boundary identification between surface and epithelial layers
- Region-of-interest selection for targeted analysis
- Physical unit conversion (pixels → micrometers)

4. Intensity Analysis
- Mean intensity calculation within ASL regions
- Intensity ratio computation between layers
- Background signal characterization

5. Batch Processing & Marking
- Multiple dataset management via tree structure
- Rectangular marker placement for analysis regions
- Selective ASL calculation within marked zones
- Comprehensive metadata tracking

Advanced Algorithms
- Epithelium Segmentation: Gradient-based boundary detection with smoothing
- Adaptive Equalization: Tile-based CLAHE with configurable parameters
- Region Filtering: Size-based and coverage-based region selection
- Curve Interpolation: Gap filling and smoothing of ASL boundaries

Installation & Setup
Prerequisites
MATLAB Requirements:
- MATLAB R2020b or later
- Image Processing Toolbox (essential)
- Signal Processing Toolbox (for CBF analysis)
- Parallel Computing Toolbox (recommended for REG tool)
