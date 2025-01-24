# Plotting utilities

arrow_head = BezierPath([
    MoveTo(Point(0, 0)),
    LineTo(Point(0.3, -1.2)),
    LineTo(Point(-0.3, -1.2)),
    ClosePath()
    ])

water_surf = RGBAf(0.0, 0.447, 0.741, 0.5)
water_bulk = RGBAf(0.0, 0.447, 0.741, 0.3)
water_bulk2 = RGBAf(0.0, 0.447, 0.741, 0.1)
sand_surf = RGBAf(0.929, 0.694, 0.125, 0.8)
sand_bulk = RGBAf(0.929, 0.694, 0.125, 0.5)
sand_bulk2 = RGBAf(0.929, 0.694, 0.125, 0.2)
plate1 = RGBAf(0.5, 0.5, 0.5, 1.0)
plate2 = RGBAf(0.5, 0.5, 0.5, 0.5)
plate3 = RGBAf(0.5, 0.5, 0.5, 0.2)
