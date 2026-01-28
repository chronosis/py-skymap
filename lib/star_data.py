"""
Static data providers for well-known bright stars, galaxies, and deep-sky objects.

These were originally defined in `skymap-gen.py` and are extracted here to
keep the main script focused on orchestration and plotting.
"""


def get_bright_galaxies():
    """Return data for the brightest distant background galaxies visible to naked eye.

    Returns list of dicts with: name, ra_deg, dec_deg, major_axis_deg, minor_axis_deg,
    position_angle_deg, apparent_mag, distance_pc.
    """
    return [
        {
            "name": "LMC",
            "ra_deg": 80.8938,  # Large Magellanic Cloud
            "dec_deg": -69.7561,
            "major_axis_deg": 10.75,  # Angular size
            "minor_axis_deg": 9.33,
            "position_angle_deg": 0,  # Approximate
            "apparent_mag": 0.9,
            "distance_pc": 50_000.0,  # ~163 kly
        },
        {
            "name": "SMC",
            "ra_deg": 13.1867,  # Small Magellanic Cloud
            "dec_deg": -72.8286,
            "major_axis_deg": 5.20,
            "minor_axis_deg": 3.25,
            "position_angle_deg": 45,  # Approximate
            "apparent_mag": 2.7,
            "distance_pc": 60_000.0,  # ~196 kly
        },
        {
            "name": "Andromeda",
            "ra_deg": 10.6847,  # M31
            "dec_deg": 41.2687,
            "major_axis_deg": 3.167,  # ~190 arcmin
            "minor_axis_deg": 1.0,  # ~60 arcmin
            "position_angle_deg": 35,  # Approximate
            "apparent_mag": 3.4,
            "distance_pc": 2_540_000.0,  # ~2.54 Mly
        },
        {
            "name": "Triangulum",
            "ra_deg": 23.4621,  # M33
            "dec_deg": 30.6602,
            "major_axis_deg": 1.0,  # ~60 arcmin
            "minor_axis_deg": 0.7,  # ~40 arcmin
            "position_angle_deg": 23,  # Approximate
            "apparent_mag": 5.7,
            "distance_pc": 3_000_000.0,  # ~3.0 Mly
        },
        {
            "name": "M3",
            "ra_deg": 205.5483,
            "dec_deg": 28.3772,
            "major_axis_deg": 0.25,
            "minor_axis_deg": 0.25,
            "position_angle_deg": 0,
            "apparent_mag": 6.2,
            "distance_pc": 10_200.0,  # ~10.2 kpc
        },
    ]


def get_bright_stars():
    """Return data for ~100 brightest and most well-known stars.

    Returns list of dicts with: name, ra_deg, dec_deg, apparent_mag, parallax_mas.
    Parallax values in milliarcseconds (mas). If None or 0, object is too distant for
    parallax measurement.
    """
    # NOTE: This list is copied verbatim from the original `skymap-gen.py`
    # so behaviour remains identical.
    return [
        # Top 20 brightest stars
        # Parallax values from Hipparcos/Gaia (mas). Distance (pc) = 1000 / parallax_mas
        {
            "name": "Sirius",
            "ra_deg": 101.2872,
            "dec_deg": -16.7161,
            "apparent_mag": -1.46,
            "parallax_mas": 379.21,
        },  # 2.64 pc
        {
            "name": "Canopus",
            "ra_deg": 95.9880,
            "dec_deg": -52.6957,
            "apparent_mag": -0.74,
            "parallax_mas": 10.43,
        },  # 95.9 pc
        {
            "name": "Rigil Kentaurus",
            "ra_deg": 219.9009,
            "dec_deg": -60.8356,
            "apparent_mag": -0.27,
            "parallax_mas": 754.81,
        },  # 1.32 pc (Alpha Centauri)
        {
            "name": "Arcturus",
            "ra_deg": 213.9153,
            "dec_deg": 19.1824,
            "apparent_mag": -0.05,
            "parallax_mas": 88.83,
        },  # 11.26 pc
        {
            "name": "Vega",
            "ra_deg": 279.2347,
            "dec_deg": 38.7837,
            "apparent_mag": 0.03,
            "parallax_mas": 130.23,
        },  # 7.68 pc
        {
            "name": "Capella",
            "ra_deg": 79.1723,
            "dec_deg": 45.9980,
            "apparent_mag": 0.08,
            "parallax_mas": 77.29,
        },  # 12.94 pc
        {
            "name": "Rigel",
            "ra_deg": 78.6345,
            "dec_deg": -8.2016,
            "apparent_mag": 0.13,
            "parallax_mas": 3.78,
        },  # 264.6 pc
        {
            "name": "Procyon",
            "ra_deg": 114.8255,
            "dec_deg": 5.2249,
            "apparent_mag": 0.34,
            "parallax_mas": 286.05,
        },  # 3.50 pc
        {
            "name": "Betelgeuse",
            "ra_deg": 88.7929,
            "dec_deg": 7.4071,
            "apparent_mag": 0.42,
            "parallax_mas": 4.51,
        },  # 221.7 pc
        {
            "name": "Achernar",
            "ra_deg": 24.4285,
            "dec_deg": -57.2368,
            "apparent_mag": 0.46,
            "parallax_mas": 23.39,
        },  # 42.75 pc
        {
            "name": "Hadar",
            "ra_deg": 210.9559,
            "dec_deg": -60.3730,
            "apparent_mag": 0.61,
            "parallax_mas": 8.32,
        },  # 120.2 pc
        {
            "name": "Altair",
            "ra_deg": 297.6958,
            "dec_deg": 8.8683,
            "apparent_mag": 0.76,
            "parallax_mas": 194.44,
        },  # 5.14 pc
        {
            "name": "Acrux",
            "ra_deg": 186.6496,
            "dec_deg": -63.0991,
            "apparent_mag": 0.77,
            "parallax_mas": 10.13,
        },  # 98.7 pc
        {
            "name": "Aldebaran",
            "ra_deg": 68.9802,
            "dec_deg": 16.5093,
            "apparent_mag": 0.85,
            "parallax_mas": 48.94,
        },  # 20.43 pc
        {
            "name": "Antares",
            "ra_deg": 247.3519,
            "dec_deg": -26.4320,
            "apparent_mag": 0.96,
            "parallax_mas": 5.89,
        },  # 169.8 pc
        {
            "name": "Spica",
            "ra_deg": 201.2983,
            "dec_deg": -11.1613,
            "apparent_mag": 0.98,
            "parallax_mas": 12.44,
        },  # 80.39 pc
        {
            "name": "Pollux",
            "ra_deg": 116.3289,
            "dec_deg": 28.0262,
            "apparent_mag": 1.14,
            "parallax_mas": 96.54,
        },  # 10.36 pc
        {
            "name": "Fomalhaut",
            "ra_deg": 344.4127,
            "dec_deg": -29.6222,
            "apparent_mag": 1.16,
            "parallax_mas": 130.08,
        },  # 7.69 pc
        {
            "name": "Deneb",
            "ra_deg": 310.3579,
            "dec_deg": 45.2803,
            "apparent_mag": 1.25,
            "parallax_mas": 2.31,
        },  # 433.0 pc
        {
            "name": "Mimosa",
            "ra_deg": 191.9303,
            "dec_deg": -59.6888,
            "apparent_mag": 1.25,
            "parallax_mas": 12.59,
        },  # 79.4 pc
        # Additional well-known stars
        {
            "name": "Polaris",
            "ra_deg": 37.9546,
            "dec_deg": 89.2641,
            "apparent_mag": 1.98,
            "parallax_mas": 7.54,
        },  # 132.6 pc
        {
            "name": "Regulus",
            "ra_deg": 152.0929,
            "dec_deg": 11.9672,
            "apparent_mag": 1.35,
            "parallax_mas": 42.09,
        },  # 23.76 pc
        {
            "name": "Castor",
            "ra_deg": 113.6494,
            "dec_deg": 31.8883,
            "apparent_mag": 1.58,
            "parallax_mas": 66.50,
        },  # 15.04 pc
        {
            "name": "Bellatrix",
            "ra_deg": 81.2828,
            "dec_deg": 6.3497,
            "apparent_mag": 1.64,
            "parallax_mas": 13.42,
        },  # 74.5 pc
        {
            "name": "Elnath",
            "ra_deg": 81.5728,
            "dec_deg": 28.6075,
            "apparent_mag": 1.65,
            "parallax_mas": 23.84,
        },  # 41.95 pc
        {
            "name": "Miaplacidus",
            "ra_deg": 138.2999,
            "dec_deg": -69.7172,
            "apparent_mag": 1.67,
            "parallax_mas": 20.71,
        },  # 48.3 pc
        {
            "name": "Alnilam",
            "ra_deg": 84.0534,
            "dec_deg": -1.2019,
            "apparent_mag": 1.69,
            "parallax_mas": 1.65,
        },  # 606 pc
        {
            "name": "Alnitak",
            "ra_deg": 85.1897,
            "dec_deg": -1.9426,
            "apparent_mag": 1.74,
            "parallax_mas": 4.43,
        },  # 225.7 pc
        {
            "name": "Mirfak",
            "ra_deg": 51.0807,
            "dec_deg": 49.8612,
            "apparent_mag": 1.79,
            "parallax_mas": 6.44,
        },  # 155.3 pc
        {
            "name": "Dubhe",
            "ra_deg": 165.932,
            "dec_deg": 61.751,
            "apparent_mag": 1.79,
            "parallax_mas": 26.54,
        },  # 37.68 pc
        {
            "name": "Wezen",
            "ra_deg": 107.0979,
            "dec_deg": -26.3932,
            "apparent_mag": 1.83,
            "parallax_mas": 2.43,
        },  # 411.5 pc
        {
            "name": "Alkaid",
            "ra_deg": 206.8852,
            "dec_deg": 49.3133,
            "apparent_mag": 1.86,
            "parallax_mas": 31.88,
        },  # 31.37 pc
        {
            "name": "Sargas",
            "ra_deg": 255.9867,
            "dec_deg": -42.9978,
            "apparent_mag": 1.87,
            "parallax_mas": 5.89,
        },  # 169.8 pc
        {
            "name": "Avior",
            "ra_deg": 125.6285,
            "dec_deg": -59.5095,
            "apparent_mag": 1.86,
            "parallax_mas": 7.51,
        },  # 133.2 pc
        {
            "name": "Menkalinan",
            "ra_deg": 89.8822,
            "dec_deg": 44.9474,
            "apparent_mag": 1.9,
            "parallax_mas": 40.16,
        },  # 24.9 pc
        {
            "name": "Atria",
            "ra_deg": 252.1662,
            "dec_deg": -69.0277,
            "apparent_mag": 1.91,
            "parallax_mas": 8.33,
        },  # 120.0 pc
        {
            "name": "Alhena",
            "ra_deg": 99.4279,
            "dec_deg": 16.5403,
            "apparent_mag": 1.93,
            "parallax_mas": 30.49,
        },  # 32.8 pc
        {
            "name": "Peacock",
            "ra_deg": 306.4119,
            "dec_deg": -56.7351,
            "apparent_mag": 1.94,
            "parallax_mas": 7.56,
        },  # 132.3 pc
        {
            "name": "Alsephina",
            "ra_deg": 140.5284,
            "dec_deg": -54.7088,
            "apparent_mag": 1.96,
            "parallax_mas": 8.46,
        },  # 118.2 pc
        {
            "name": "Mirzam",
            "ra_deg": 95.6749,
            "dec_deg": -17.9559,
            "apparent_mag": 1.98,
            "parallax_mas": 3.29,
        },  # 304.0 pc
        {
            "name": "Alphard",
            "ra_deg": 141.8968,
            "dec_deg": -8.6586,
            "apparent_mag": 1.99,
            "parallax_mas": 41.35,
        },  # 24.18 pc
        {
            "name": "Algieba",
            "ra_deg": 154.9926,
            "dec_deg": 19.8415,
            "apparent_mag": 2.01,
            "parallax_mas": 25.96,
        },  # 38.52 pc
        {
            "name": "Diphda",
            "ra_deg": 10.8974,
            "dec_deg": -17.9866,
            "apparent_mag": 2.04,
            "parallax_mas": 33.62,
        },  # 29.74 pc
        {
            "name": "Mizar",
            "ra_deg": 200.9814,
            "dec_deg": 54.9254,
            "apparent_mag": 2.04,
            "parallax_mas": 41.73,
        },  # 23.96 pc
        {
            "name": "Nunki",
            "ra_deg": 283.8164,
            "dec_deg": -26.2961,
            "apparent_mag": 2.05,
            "parallax_mas": 13.87,
        },  # 72.1 pc
        {
            "name": "Kaus Australis",
            "ra_deg": 276.043,
            "dec_deg": -34.3846,
            "apparent_mag": 1.79,
            "parallax_mas": 8.34,
        },  # 119.9 pc
        {
            "name": "Sadr",
            "ra_deg": 305.5571,
            "dec_deg": 40.2567,
            "apparent_mag": 2.23,
            "parallax_mas": 1.78,
        },  # 561.8 pc
        {
            "name": "Eltanin",
            "ra_deg": 262.6082,
            "dec_deg": 51.4889,
            "apparent_mag": 2.24,
            "parallax_mas": 21.09,
        },  # 47.4 pc
        {
            "name": "Kaus Media",
            "ra_deg": 274.4067,
            "dec_deg": -29.8281,
            "apparent_mag": 2.7,
            "parallax_mas": 8.59,
        },  # 116.4 pc
        {
            "name": "Alpheratz",
            "ra_deg": 2.0969,
            "dec_deg": 29.0904,
            "apparent_mag": 2.07,
            "parallax_mas": 33.62,
        },  # 29.74 pc
        {
            "name": "Mirach",
            "ra_deg": 17.433,
            "dec_deg": 35.6206,
            "apparent_mag": 2.07,
            "parallax_mas": 16.52,
        },  # 60.5 pc
        {
            "name": "Rasalgethi",
            "ra_deg": 258.6619,
            "dec_deg": 14.3903,
            "apparent_mag": 2.78,
            "parallax_mas": 8.64,
        },  # 115.7 pc
        {
            "name": "Kochab",
            "ra_deg": 222.6764,
            "dec_deg": 74.1555,
            "apparent_mag": 2.08,
            "parallax_mas": 17.49,
        },  # 57.2 pc
        {
            "name": "Saiph",
            "ra_deg": 86.9391,
            "dec_deg": -9.6696,
            "apparent_mag": 2.07,
            "parallax_mas": 5.04,
        },  # 198.4 pc
        {
            "name": "Hamal",
            "ra_deg": 31.7934,
            "dec_deg": 23.4624,
            "apparent_mag": 2.01,
            "parallax_mas": 49.56,
        },  # 20.18 pc
        {
            "name": "Algol",
            "ra_deg": 47.0422,
            "dec_deg": 40.9556,
            "apparent_mag": 2.09,
            "parallax_mas": 35.14,
        },  # 28.46 pc
        {
            "name": "Dschubba",
            "ra_deg": 240.0833,
            "dec_deg": -22.6217,
            "apparent_mag": 2.29,
            "parallax_mas": 5.48,
        },  # 182.5 pc
        {
            "name": "Zubeneschamali",
            "ra_deg": 229.2517,
            "dec_deg": -9.3829,
            "apparent_mag": 2.61,
            "parallax_mas": 18.77,
        },  # 53.3 pc
        {
            "name": "Graffias",
            "ra_deg": 244.5804,
            "dec_deg": -26.4321,
            "apparent_mag": 2.62,
            "parallax_mas": 5.89,
        },  # 169.8 pc
        {
            "name": "Iota Carinae",
            "ra_deg": 139.2725,
            "dec_deg": -59.2752,
            "apparent_mag": 2.21,
            "parallax_mas": 4.2,
        },  # 238.1 pc
        {
            "name": "Theta Carinae",
            "ra_deg": 160.7392,
            "dec_deg": -64.3944,
            "apparent_mag": 2.74,
            "parallax_mas": 5.55,
        },  # 180.2 pc
        {
            "name": "Aspidiske",
            "ra_deg": 141.5269,
            "dec_deg": -64.3944,
            "apparent_mag": 2.21,
            "parallax_mas": 5.4,
        },  # 185.2 pc
        # Additional bright stars to reach ~100
        {
            "name": "Schedar",
            "ra_deg": 10.1268,
            "dec_deg": 56.5373,
            "apparent_mag": 2.24,
            "parallax_mas": 14.29,
        },  # 70.0 pc
        {
            "name": "Caph",
            "ra_deg": 2.2945,
            "dec_deg": 59.1498,
            "apparent_mag": 2.28,
            "parallax_mas": 60.42,
        },  # 16.55 pc
        {
            "name": "Achird",
            "ra_deg": 9.8322,
            "dec_deg": 54.2844,
            "apparent_mag": 3.46,
            "parallax_mas": 168.45,
        },  # 5.94 pc
        {
            "name": "Almach",
            "ra_deg": 30.9748,
            "dec_deg": 42.3297,
            "apparent_mag": 2.1,
            "parallax_mas": 16.64,
        },  # 60.1 pc
        {
            "name": "Mira",
            "ra_deg": 34.8367,
            "dec_deg": -2.9774,
            "apparent_mag": 2.0,
            "parallax_mas": 10.91,
        },  # 91.7 pc
        {
            "name": "Menkar",
            "ra_deg": 45.5699,
            "dec_deg": 4.0897,
            "apparent_mag": 2.54,
            "parallax_mas": 15.79,
        },  # 63.3 pc
        {
            "name": "Baten Kaitos",
            "ra_deg": 13.6605,
            "dec_deg": -10.335,
            "apparent_mag": 3.74,
            "parallax_mas": 16.23,
        },  # 61.6 pc
        {
            "name": "Ankaa",
            "ra_deg": 6.5708,
            "dec_deg": -42.3058,
            "apparent_mag": 2.4,
            "parallax_mas": 40.9,
        },  # 24.45 pc
        {
            "name": "Markab",
            "ra_deg": 346.1902,
            "dec_deg": 15.2053,
            "apparent_mag": 2.49,
            "parallax_mas": 15.07,
        },  # 66.4 pc
        {
            "name": "Scheat",
            "ra_deg": 345.9436,
            "dec_deg": 28.0828,
            "apparent_mag": 2.44,
            "parallax_mas": 16.64,
        },  # 60.1 pc
        {
            "name": "Algenib",
            "ra_deg": 3.3089,
            "dec_deg": 15.1836,
            "apparent_mag": 2.83,
            "parallax_mas": 12.1,
        },  # 82.6 pc
        {
            "name": "Enif",
            "ra_deg": 326.0465,
            "dec_deg": 9.875,
            "apparent_mag": 2.38,
            "parallax_mas": 4.97,
        },  # 201.2 pc
        {
            "name": "Homam",
            "ra_deg": 340.7508,
            "dec_deg": 10.8314,
            "apparent_mag": 3.4,
            "parallax_mas": 20.23,
        },  # 49.4 pc
        {
            "name": "Matar",
            "ra_deg": 330.68,
            "dec_deg": 30.2211,
            "apparent_mag": 2.99,
            "parallax_mas": 13.33,
        },  # 75.0 pc
        {
            "name": "Biham",
            "ra_deg": 340.3652,
            "dec_deg": 6.1978,
            "apparent_mag": 3.51,
            "parallax_mas": 20.23,
        },  # 49.4 pc
        {
            "name": "Sadalbari",
            "ra_deg": 344.4127,
            "dec_deg": -29.6222,
            "apparent_mag": 3.51,
            "parallax_mas": 22.46,
        },  # 44.5 pc
    ]


def get_bright_deep_sky_objects():
    """Return data for magnitude 7 or brighter nebulas, globular clusters, and other non-star objects.

    Returns list of dicts with: name, ra_deg, dec_deg, major_axis_deg, minor_axis_deg,
    position_angle_deg, apparent_mag, object_type, distance_pc.
    """
    # NOTE: Copied directly from original script; no behavioural changes.
    return [
        # Bright Nebulas
        {
            "name": "Orion Nebula",
            "ra_deg": 83.8221,
            "dec_deg": -5.3911,
            "major_axis_deg": 1.5,
            "minor_axis_deg": 1.0,
            "position_angle_deg": 0,
            "apparent_mag": 4.0,
            "object_type": "nebula",
            "distance_pc": 414.0,
        },  # VLBA parallax
        {
            "name": "Carina Nebula",
            "ra_deg": 160.8950,
            "dec_deg": -59.6856,
            "major_axis_deg": 2.0,
            "minor_axis_deg": 2.0,
            "position_angle_deg": 0,
            "apparent_mag": 1.0,
            "object_type": "nebula",
            "distance_pc": 2350.0,
        },  # Gaia EDR3
        {
            "name": "Eagle Nebula",
            "ra_deg": 274.7,
            "dec_deg": -13.8067,
            "major_axis_deg": 0.5,
            "minor_axis_deg": 0.5,
            "position_angle_deg": 0,
            "apparent_mag": 6.0,
            "object_type": "nebula",
            "distance_pc": 2000.0,
        },  # ~2 kpc
        {
            "name": "Lagoon Nebula",
            "ra_deg": 271.0,
            "dec_deg": -24.3833,
            "major_axis_deg": 1.33,
            "minor_axis_deg": 1.0,
            "position_angle_deg": 0,
            "apparent_mag": 6.0,
            "object_type": "nebula",
            "distance_pc": 4100.0,
        },  # ~4.1 kpc
        {
            "name": "Trifid Nebula",
            "ra_deg": 270.67,
            "dec_deg": -23.0167,
            "major_axis_deg": 0.5,
            "minor_axis_deg": 0.5,
            "position_angle_deg": 0,
            "apparent_mag": 6.3,
            "object_type": "nebula",
            "distance_pc": 4100.0,
        },  # ~4.1 kpc (near Lagoon)
        {
            "name": "Omega Nebula",
            "ra_deg": 275.2,
            "dec_deg": -16.15,
            "major_axis_deg": 0.67,
            "minor_axis_deg": 0.5,
            "position_angle_deg": 0,
            "apparent_mag": 6.0,
            "object_type": "nebula",
            "distance_pc": 2000.0,
        },  # ~2 kpc
        {
            "name": "Rosette Nebula",
            "ra_deg": 97.95,
            "dec_deg": 4.95,
            "major_axis_deg": 1.33,
            "minor_axis_deg": 1.0,
            "position_angle_deg": 0,
            "apparent_mag": 4.8,
            "object_type": "nebula",
            "distance_pc": 1600.0,
        },  # ~1.6 kpc
        {
            "name": "Horsehead Nebula",
            "ra_deg": 85.25,
            "dec_deg": -2.45,
            "major_axis_deg": 0.17,
            "minor_axis_deg": 0.17,
            "position_angle_deg": 0,
            "apparent_mag": 6.8,
            "object_type": "nebula",
            "distance_pc": 414.0,
        },  # Same region as Orion
        {
            "name": "North America Nebula",
            "ra_deg": 314.7,
            "dec_deg": 44.0167,
            "major_axis_deg": 2.0,
            "minor_axis_deg": 1.5,
            "position_angle_deg": 0,
            "apparent_mag": 4.0,
            "object_type": "nebula",
            "distance_pc": 1800.0,
        },  # ~1.8 kpc
        {
            "name": "Veil Nebula",
            "ra_deg": 312.75,
            "dec_deg": 30.7833,
            "major_axis_deg": 3.0,
            "minor_axis_deg": 2.5,
            "position_angle_deg": 0,
            "apparent_mag": 7.0,
            "object_type": "nebula",
            "distance_pc": 1470.0,
        },  # ~1.47 kpc (supernova remnant)
        {
            "name": "Dumbbell Nebula",
            "ra_deg": 299.9017,
            "dec_deg": 22.7211,
            "major_axis_deg": 0.25,
            "minor_axis_deg": 0.17,
            "position_angle_deg": 0,
            "apparent_mag": 7.4,
            "object_type": "nebula",
            "distance_pc": 1360.0,
        },  # Planetary nebula
        {
            "name": "Ring Nebula",
            "ra_deg": 283.3961,
            "dec_deg": 33.0292,
            "major_axis_deg": 0.08,
            "minor_axis_deg": 0.08,
            "position_angle_deg": 0,
            "apparent_mag": 8.8,
            "object_type": "nebula",
            "distance_pc": 2300.0,
        },  # Planetary nebula
        {
            "name": "Helix Nebula",
            "ra_deg": 337.41,
            "dec_deg": -20.8333,
            "major_axis_deg": 0.5,
            "minor_axis_deg": 0.5,
            "position_angle_deg": 0,
            "apparent_mag": 7.3,
            "object_type": "nebula",
            "distance_pc": 695.0,
        },  # ~695 pc (planetary nebula)
        {
            "name": "Crab Nebula",
            "ra_deg": 83.6331,
            "dec_deg": 22.0144,
            "major_axis_deg": 0.17,
            "minor_axis_deg": 0.17,
            "position_angle_deg": 0,
            "apparent_mag": 8.4,
            "object_type": "nebula",
            "distance_pc": 2000.0,
        },  # ~2 kpc (supernova remnant)
        # Globular Clusters
        {
            "name": "Omega Centauri",
            "ra_deg": 201.6967,
            "dec_deg": -47.4794,
            "major_axis_deg": 0.55,
            "minor_axis_deg": 0.55,
            "position_angle_deg": 0,
            "apparent_mag": 3.7,
            "object_type": "globular",
            "distance_pc": 5200.0,
        },  # ~5.2 kpc
        {
            "name": "47 Tucanae",
            "ra_deg": 6.0229,
            "dec_deg": -72.0814,
            "major_axis_deg": 0.5,
            "minor_axis_deg": 0.5,
            "position_angle_deg": 0,
            "apparent_mag": 4.0,
            "object_type": "globular",
            "distance_pc": 4500.0,
        },  # ~4.5 kpc
        {
            "name": "M13",
            "ra_deg": 250.4233,
            "dec_deg": 36.4614,
            "major_axis_deg": 0.33,
            "minor_axis_deg": 0.33,
            "position_angle_deg": 0,
            "apparent_mag": 5.8,
            "object_type": "globular",
            "distance_pc": 7100.0,
        },  # ~7.1 kpc
        {
            "name": "M3",
            "ra_deg": 205.5483,
            "dec_deg": 28.3772,
            "major_axis_deg": 0.25,
            "minor_axis_deg": 0.25,
            "position_angle_deg": 0,
            "apparent_mag": 6.2,
            "object_type": "globular",
            "distance_pc": 10_200.0,
        },  # ~10.2 kpc
        {
            "name": "M5",
            "ra_deg": 229.6383,
            "dec_deg": 2.0811,
            "major_axis_deg": 0.25,
            "minor_axis_deg": 0.25,
            "position_angle_deg": 0,
            "apparent_mag": 5.6,
            "object_type": "globular",
            "distance_pc": 7500.0,
        },  # ~7.5 kpc
        {
            "name": "M4",
            "ra_deg": 245.8967,
            "dec_deg": -26.5256,
            "major_axis_deg": 0.33,
            "minor_axis_deg": 0.33,
            "position_angle_deg": 0,
            "apparent_mag": 5.6,
            "object_type": "globular",
            "distance_pc": 2200.0,
        },  # ~2.2 kpc (closest globular)
        {
            "name": "M22",
            "ra_deg": 279.1,
            "dec_deg": -23.9047,
            "major_axis_deg": 0.33,
            "minor_axis_deg": 0.33,
            "position_angle_deg": 0,
            "apparent_mag": 5.1,
            "object_type": "globular",
            "distance_pc": 3200.0,
        },  # ~3.2 kpc
        {
            "name": "M15",
            "ra_deg": 322.4933,
            "dec_deg": 12.1672,
            "major_axis_deg": 0.25,
            "minor_axis_deg": 0.25,
            "position_angle_deg": 0,
            "apparent_mag": 6.2,
            "object_type": "globular",
            "distance_pc": 10_400.0,
        },  # ~10.4 kpc
        {
            "name": "M2",
            "ra_deg": 323.3625,
            "dec_deg": -0.8233,
            "major_axis_deg": 0.17,
            "minor_axis_deg": 0.17,
            "position_angle_deg": 0,
            "apparent_mag": 6.5,
            "object_type": "globular",
            "distance_pc": 11_500.0,
        },  # ~11.5 kpc
        {
            "name": "M92",
            "ra_deg": 259.2817,
            "dec_deg": 43.1358,
            "major_axis_deg": 0.17,
            "minor_axis_deg": 0.17,
            "position_angle_deg": 0,
            "apparent_mag": 6.4,
            "object_type": "globular",
            "distance_pc": 8300.0,
        },  # ~8.3 kpc
        # Open Clusters
        {
            "name": "Pleiades",
            "ra_deg": 56.8711,
            "dec_deg": 24.1053,
            "major_axis_deg": 2.0,
            "minor_axis_deg": 2.0,
            "position_angle_deg": 0,
            "apparent_mag": 1.6,
            "object_type": "cluster",
            "distance_pc": 135.0,
        },  # Hipparcos/Gaia parallax
        {
            "name": "Hyades",
            "ra_deg": 66.7325,
            "dec_deg": 15.87,
            "major_axis_deg": 5.5,
            "minor_axis_deg": 4.0,
            "position_angle_deg": 0,
            "apparent_mag": 0.5,
            "object_type": "cluster",
            "distance_pc": 46.35,
        },  # Hipparcos parallax
        {
            "name": "Beehive Cluster",
            "ra_deg": 130.1,
            "dec_deg": 19.6667,
            "major_axis_deg": 1.5,
            "minor_axis_deg": 1.5,
            "position_angle_deg": 0,
            "apparent_mag": 3.7,
            "object_type": "cluster",
            "distance_pc": 187.0,
        },  # M44, ~187 pc
        {
            "name": "Double Cluster",
            "ra_deg": 34.75,
            "dec_deg": 57.15,
            "major_axis_deg": 1.0,
            "minor_axis_deg": 1.0,
            "position_angle_deg": 0,
            "apparent_mag": 4.3,
            "object_type": "cluster",
            "distance_pc": 2300.0,
        },  # h+Chi Persei, ~2.3 kpc
        {
            "name": "M6",
            "ra_deg": 265.0833,
            "dec_deg": -32.2167,
            "major_axis_deg": 1.0,
            "minor_axis_deg": 1.0,
            "position_angle_deg": 0,
            "apparent_mag": 4.2,
            "object_type": "cluster",
            "distance_pc": 1600.0,
        },  # ~1.6 kpc
        {
            "name": "M7",
            "ra_deg": 268.4708,
            "dec_deg": -34.7917,
            "major_axis_deg": 1.33,
            "minor_axis_deg": 1.33,
            "position_angle_deg": 0,
            "apparent_mag": 3.3,
            "object_type": "cluster",
            "distance_pc": 980.0,
        },  # ~980 pc
        {
            "name": "M11",
            "ra_deg": 282.775,
            "dec_deg": -6.2667,
            "major_axis_deg": 0.67,
            "minor_axis_deg": 0.67,
            "position_angle_deg": 0,
            "apparent_mag": 5.8,
            "object_type": "cluster",
            "distance_pc": 6100.0,
        },  # ~6.1 kpc
        {
            "name": "M44",
            "ra_deg": 130.1,
            "dec_deg": 19.6667,
            "major_axis_deg": 1.5,
            "minor_axis_deg": 1.5,
            "position_angle_deg": 0,
            "apparent_mag": 3.7,
            "object_type": "cluster",
            "distance_pc": 187.0,
        },  # Beehive Cluster, same as above
        {
            "name": "M45",
            "ra_deg": 56.8711,
            "dec_deg": 24.1053,
            "major_axis_deg": 2.0,
            "minor_axis_deg": 2.0,
            "position_angle_deg": 0,
            "apparent_mag": 1.6,
            "object_type": "cluster",
            "distance_pc": 135.0,
        },  # Pleiades, same as above
    ]

