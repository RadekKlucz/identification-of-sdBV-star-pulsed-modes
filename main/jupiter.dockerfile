# The base image 
FROM jupyter/base-notebook

# Working directory
WORKDIR /home/notebooks

# Requirements file
COPY --chown=${NB_UID}:${NB_GID} requirements.txt . 

# Install requirements module
RUN pip install --quiet --no-cache-dir --requirement requirements.txt && \
    fix-permissions "${CONDA_DIR}" && fix-permissions "/home/${NB_USER}"

# Environmental variable
ENV JUPYTER_ENABLE_LAB=yes

# Run shell command for notebook on start 
CMD jupyter notebook --port=8888 --no-browser --ip=0.0.0.0 --allow-root