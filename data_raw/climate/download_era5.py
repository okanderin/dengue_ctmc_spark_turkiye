import cdsapi

client = cdsapi.Client()

client.retrieve(
    'reanalysis-era5-land-monthly-means',
    {
        'product_type': 'monthly_averaged_reanalysis',
        'variable': [
            '2m_temperature',
            '2m_dewpoint_temperature',
            'total_precipitation',
        ],
        'year': [str(y) for y in range(1981, 2015)],
        'month': [f'{m:02d}' for m in range(1, 13)],
        'time': '00:00',
        'area': [42.5, 25.5, 35.5, 45.5],
        'format': 'netcdf',
    },
    'era5_land_1981_2014.nc'
)
