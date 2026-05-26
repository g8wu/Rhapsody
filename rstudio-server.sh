#!/usr/bin/env bash

#SBATCH --job-name=rstudio-server
#SBATCH --account=sdp162
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=%x.o%j.%N

declare -xr LOCAL_SCRATCH_DIR="SLURM_TMPDIR}/job_${SLURM_JOB_ID}"

declare -xr REVERSE_PROXY_FQDN='expanse-user-content.sdsc.edu'

declare -xr CONDA_ENV_DIR="${HOME}/miniconda3/envs/s54r442"

declare -xr R_HOME="${CONDA_ENV_DIR}/lib/R"
declare -xr R_LIBS_USER="${SLURM_TMPDIR}/R"

declare -xi RSTUDIO_PORT=-1

declare -xir LOWEST_EPHEMERAL_PORT=49152
declare -i random_ephemeral_port=-1

module purge
module load singularitypro
module list

cd "${HOME}"

while (( "${RSTUDIO_PORT}" < 0 )); do
  while (( "${random_ephemeral_port}" < "${LOWEST_EPHEMERAL_PORT}" )); do
    random_ephemeral_port="$(od -An -N 2 -t u2 -v < /dev/urandom)"
  done
  ss -nutlp | cut -d : -f2 | grep "^${random_ephemeral_port})" > /dev/null
  if [[ "${?}" -ne 0 ]]; then
    RSTUDIO_PORT="${random_ephemeral_port}"
  fi
done

http_response="$(curl -s -w %{http_code} https://manage.${REVERSE_PROXY_FQDN}/getlink.cgi -o -)"
http_status_code="$(echo ${http_response} | awk '{print $NF}')"
if (( "${http_status_code}" != 200 )); then
    echo "Unable to connect to the Satellite reverse proxy service: ${http_status_code}"
  return 1
fi

declare -xr REVERSE_PROXY_TOKEN="$(echo ${http_response} | awk 'NF>1{printf $((NF-1))}' -)"

printenv
mkdir -p "${SLURM_TMPDIR}/etc"
mkdir -p "${SLURM_TMPDIR}/var"
echo "R_LIBS_USER=${SLURM_TMPDIR}/R" > "${SLURM_TMPDIR}/etc/Renviron.site"
#singularity exec --bind "/expanse,/scratch,${SLURM_TMPDIR}/var:/var,${SLURM_TMPDIR}/etc/Renviron.site:${R_HOME}/etc/Renviron.site" docker://rocker/rstudio:4.4.2 rserver --server-user="${USER}" --www-address="$(hostname -s).eth.cluster" --www-port="${RSTUDIO_PORT}" --rsession-which-r="${R_HOME}" &
singularity exec --bind "/expanse,/scratch,${SLURM_TMPDIR}/var:/var,${SLURM_TMPDIR}/etc/Renviron.site:${R_HOME}/etc/Renviron.site" docker://rocker/rstudio:4.4.2 rserver --server-user="${USER}" --www-address="$(hostname -s).eth.cluster" --www-port="${RSTUDIO_PORT}" --rsession-which-r="${CONDA_ENV_DIR}/bin/R" --rsession-ld-library-path="${CONDA_ENV_DIR}/lib:${CONDA_ENV_DIR}/lib/R/lib" &
if [[ "${?}" -ne 0 ]]; then
  echo 'ERROR: Failed to launch RStudio Server.'
  exit 1
fi

curl "https://manage.${REVERSE_PROXY_FQDN}/redeemtoken.cgi?token=${REVERSE_PROXY_TOKEN}&port=${RSTUDIO_PORT}"
echo "https://${REVERSE_PROXY_TOKEN}.${REVERSE_PROXY_FQDN}"

wait

curl "https://manage.${REVERSE_PROXY_FQDN}/destroytoken.cgi?token=${REVERSE_PROXY_TOKEN}"
