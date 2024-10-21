# Simulated Data

This folder contains simulated data for the ray tracer 3D project. The data is organized by distance from the detector and thickness of the oil.

## Folder Structure

- **Distance Folder**: Each folder represents a different distance from the detector.
    - **Thickness Folder**: Within each distance folder, there are subfolders for each thickness of the oil.
        - `data.txt`: Contains the raw data.
        - `plot.png`: Contains the plot of the data.
    - `total_plot.png`: Contains the total plot for all thicknesses at the given distance.

## Example Structure

```
/simulated_data/
    /distance_1/
        /thickness_1/
            data.txt
            plot.png
        /thickness_2/
            data.txt
            plot.png
        total_plot.png
    /distance_2/
        /thickness_1/
            data.txt
            plot.png
        /thickness_2/
            data.txt
            plot.png
        total_plot.png
```

## Usage

Use the provided data and plots to analyze the behavior of the ray tracer at different distances and oil thicknesses.
