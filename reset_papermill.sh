killall fish_queue_processor.sh
killall -u simone python
ssh monod10 'rm -Rf /tmp/dask-worker-space && killall -u simone python'
ssh monod11 'rm -Rf /tmp/dask-worker-space && killall -u simone python'
ssh monod12 'rm -Rf /tmp/dask-worker-space && killall -u simone python'
ssh monod33 'rm -Rf /tmp/dask-worker-space && killall -u simone python'