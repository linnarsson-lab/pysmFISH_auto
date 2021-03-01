
def nd2_parsing 

if parsing_type == 'original'
    # ----------------------------------------------------------------
    # PARSING THE MICROSCOPY DATA
    start = time.time()
    logger.info(f'start reparsing raw data')
    # Create empty zarr file for the parse data
    parsed_raw_data_fpath = create_empty_zarr_file(experiment_fpath=experiment_fpath,
                                        tag=parsed_image_tag)

    # Parse the data
    all_raw_nd2 = nd2_raw_files_selector(experiment_fpath)

    parsing_futures = client.map(nikon_nd2_autoparser_zarr,
                                all_raw_nd2,
                                parsed_raw_data_fpath=parsed_raw_data_fpath,
                                experiment_info=experiment_info)

    _ = client.gather(parsing_futures)

    logger.info(f'reparsing completed in {(time.time()-start)/60} min')
    # ----------------------------------------------------------------