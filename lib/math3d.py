import numpy as np
from astropy.coordinates import SkyCoord, Distance, CartesianRepresentation
from astropy import units as u


def get_relative_coords(target_xyz, object_xyz):
    """Compute vector, distance (pc), and unit vector from target to object.

    Handles both scalar (3,) and array (N, 3) inputs for object_xyz.
    Returns:
        vector_from_target: object_xyz - target_xyz
        distance_pc: Euclidean distance in parsecs
        unit_vector: normalized direction vectors
    """
    target_xyz = np.asarray(target_xyz)
    object_xyz = np.asarray(object_xyz)

    # Ensure proper shapes for broadcasting
    if target_xyz.ndim == 1 and object_xyz.ndim == 2:
        # (3,) and (N,3) broadcast naturally
        pass
    elif target_xyz.ndim == 2 and object_xyz.ndim == 1:
        object_xyz = object_xyz.reshape(1, -1)

    vector_from_target = object_xyz - target_xyz
    distance = np.linalg.norm(vector_from_target, axis=-1, keepdims=False)

    # Avoid division by zero
    safe_distance = np.where(distance > 1e-10, distance, 1.0)
    unit_vector = vector_from_target / np.expand_dims(safe_distance, -1)

    return vector_from_target, distance, unit_vector


def transform_to_target_frame(
    target_coord: SkyCoord,
    object_ra_deg: float,
    object_dec_deg: float,
    object_distance_pc: float,
    use_earth_centric_approx: bool = False,
):
    """Transform an object's coordinates into the target star's local frame.

    Returns (azimuth_rad, elevation_rad, z_component).
    """
    if use_earth_centric_approx:
        object_icrs = SkyCoord(
            ra=object_ra_deg * u.deg,
            dec=object_dec_deg * u.deg,
            distance=Distance(parallax=1.0 * u.mas),
            frame="icrs",
        )
    else:
        object_icrs = SkyCoord(
            ra=object_ra_deg * u.deg,
            dec=object_dec_deg * u.deg,
            distance=object_distance_pc * u.pc,
            frame="icrs",
        )

    object_cart = object_icrs.cartesian.xyz.value
    target_cart = target_coord.cartesian.xyz.value

    _, _, unit_vector = get_relative_coords(target_cart, object_cart)

    x, y, z = unit_vector
    azimuth_rad = np.arctan2(y, x)
    azimuth_rad = np.mod(azimuth_rad, 2 * np.pi)
    elevation_rad = np.arcsin(np.clip(z, -1.0, 1.0))

    return azimuth_rad, elevation_rad, z

