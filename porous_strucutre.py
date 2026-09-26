import math
from pathlib import Path

import cadquery as cq
from cadquery import exporters


# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

SAMPLE_DIAMETER = 28.0   # mm
SAMPLE_HEIGHT = 54.0     # mm, flow direction length
CHANNEL_SPACING = 5.0    # mm

TARGET_POROSITIES = [0.01, 0.025, 0.05]

OUTPUT_DIR = "porous_cylinders"

# STL export quality
STL_TOLERANCE = 0.02          # mm, smaller = finer STL
STL_ANGULAR_TOLERANCE = 0.1   # rad, smaller = finer round surfaces


# ------------------------------------------------------------
# Geometry helper functions
# ------------------------------------------------------------

def get_channel_centers(
    sample_diameter: float,
    channel_spacing: float,
    max_channel_radius: float,
    edge_clearance: float = 0.2,
):
    """
    Generate channel center positions inside the circular cross-section.

    The centers are chosen so that even the largest channel radius still
    remains inside the sample boundary.
    """

    sample_radius = sample_diameter / 2.0

    n = int(math.floor(sample_radius / channel_spacing))
    positions = [i * channel_spacing for i in range(-n, n + 1)]

    centers = []

    for x in positions:
        for y in positions:
            distance_from_center = math.sqrt(x**2 + y**2)

            if distance_from_center + max_channel_radius + edge_clearance <= sample_radius:
                centers.append((x, y))

    if len(centers) == 0:
        raise ValueError("No channel centers found. Reduce spacing or channel radius.")

    return centers


def calculate_radius_from_porosity(
    target_porosity: float,
    sample_diameter: float,
    number_of_channels: int,
):
    """
    For straight, parallel, non-overlapping axial channels:

        phi = N * r^2 / R^2

    Therefore:

        r = R * sqrt(phi / N)
    """

    sample_radius = sample_diameter / 2.0

    channel_radius = sample_radius * math.sqrt(
        target_porosity / number_of_channels
    )

    return channel_radius


def create_channel_solid(
    channel_centers,
    channel_radius: float,
    sample_height: float,
):
    """
    Create axial cylindrical channels in z-direction.
    These are the negative volumes that are later subtracted.
    """

    channel_length = sample_height * 1.5
    channel_solids = []

    for x, y in channel_centers:
        cyl = cq.Solid.makeCylinder(
            channel_radius,
            channel_length,
            cq.Vector(x, y, -channel_length / 2.0),
            cq.Vector(0, 0, 1),
        )
        channel_solids.append(cyl)

    return cq.Compound.makeCompound(channel_solids)


def create_porous_cylinder(
    sample_diameter: float,
    sample_height: float,
    channel_radius: float,
    channel_centers,
):
    """
    Create solid cylinder and subtract axial channels.
    """

    sample_radius = sample_diameter / 2.0

    cylinder = (
        cq.Workplane("XY")
        .cylinder(sample_height, sample_radius)
        .val()
    )

    channels = create_channel_solid(
        channel_centers=channel_centers,
        channel_radius=channel_radius,
        sample_height=sample_height,
    )

    porous_cylinder = cylinder.cut(channels)

    return porous_cylinder


def calculate_porosity_cylinder(
    shape,
    sample_diameter: float,
    sample_height: float,
):
    """
    Calculate CAD porosity from the remaining solid volume.

        phi = 1 - V_solid / V_total
    """

    sample_radius = sample_diameter / 2.0

    total_volume = math.pi * sample_radius**2 * sample_height
    solid_volume = shape.Volume()

    porosity = 1.0 - solid_volume / total_volume

    return porosity


# ------------------------------------------------------------
# Main script
# ------------------------------------------------------------

def main():
    output_path = Path(OUTPUT_DIR)
    output_path.mkdir(exist_ok=True)

    # Choose a maximum radius for defining stable channel positions.
    # 0.45 * spacing keeps a small wall between neighboring channels.
    max_channel_radius = 0.45 * CHANNEL_SPACING

    channel_centers = get_channel_centers(
        sample_diameter=SAMPLE_DIAMETER,
        channel_spacing=CHANNEL_SPACING,
        max_channel_radius=max_channel_radius,
        edge_clearance=0.2,
    )

    number_of_channels = len(channel_centers)

    print("----------------------------------------")
    print("Porous cylinder generation")
    print("----------------------------------------")
    print(f"Sample diameter:     {SAMPLE_DIAMETER:.2f} mm")
    print(f"Sample height:       {SAMPLE_HEIGHT:.2f} mm")
    print(f"Channel spacing:     {CHANNEL_SPACING:.2f} mm")
    print(f"Number of channels:  {number_of_channels}")
    print(f"Max channel radius:  {max_channel_radius:.3f} mm")
    print("----------------------------------------")

    for target_phi in TARGET_POROSITIES:
        channel_radius = calculate_radius_from_porosity(
            target_porosity=target_phi,
            sample_diameter=SAMPLE_DIAMETER,
            number_of_channels=number_of_channels,
        )

        if channel_radius > max_channel_radius:
            print()
            print(f"Target porosity {target_phi:.2f} cannot be reached.")
            print(f"Required radius: {channel_radius:.3f} mm")
            print(f"Maximum radius:  {max_channel_radius:.3f} mm")
            print("Use smaller channel spacing or allow thinner walls.")
            continue

        shape = create_porous_cylinder(
            sample_diameter=SAMPLE_DIAMETER,
            sample_height=SAMPLE_HEIGHT,
            channel_radius=channel_radius,
            channel_centers=channel_centers,
        )

        cad_phi = calculate_porosity_cylinder(
            shape=shape,
            sample_diameter=SAMPLE_DIAMETER,
            sample_height=SAMPLE_HEIGHT,
        )

        filename = (
            f"porous_cylinder_"
            f"D{SAMPLE_DIAMETER:.0f}mm_"
            f"H{SAMPLE_HEIGHT:.0f}mm_"
            f"phi{cad_phi:.3f}_"
            f"d{2 * channel_radius:.3f}mm.stl"
        )

        filepath = output_path / filename

        exporters.export(
            cq.Workplane("XY").add(shape),
            str(filepath),
            tolerance=STL_TOLERANCE,
            angularTolerance=STL_ANGULAR_TOLERANCE,
        )

        print()
        print(f"Target porosity:     {target_phi:.3f}")
        print(f"CAD porosity:        {cad_phi:.3f}")
        print(f"Channel radius:      {channel_radius:.3f} mm")
        print(f"Channel diameter:    {2 * channel_radius:.3f} mm")
        print(f"Exported STL:        {filepath}")


if __name__ == "__main__":
    main()