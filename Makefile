APP = imicrobe-16SrDNA-OTU-Clustering
VERSION = 0.0.1
EMAIL = jklynch@email.arizona.edu

clean:
	find . \( -name \*.conf -o -name \*.out -o -name \*.log -o -name \*.param -o -name launcher_jobfile_\* \) -exec rm {} \;

files-delete:
	files-delete $(CYVERSEUSERNAME)/applications/$(APP)-$(VERSION)

files-upload:
	files-upload -F stampede/ $(CYVERSEUSERNAME)/applications/$(APP)-$(VERSION)

apps-addupdate:
	apps-addupdate -F stampede/app.json

deploy-app: clean files-delete files-upload apps-addupdate

share-app:
	systems-roles-addupdate -v -u $(TARGETUSER) -r USER tacc-stampede-$(CYVERSEUSERNAME)
	apps-pems-update -v -u $(TARGETUSER) -p READ_EXECUTE $(APP)-$(VERSION)

test:
	sbatch test.sh

jobs-submit:
	jobs-submit -F stampede/job.json

container:
	rm -f singularity/$(APP).img
	sudo singularity create --size 2048 singularity/$(APP).img
	sudo singularity bootstrap singularity/$(APP).img singularity/Singularity
	sudo chown --reference=singularity/Singularity singularity/${APP}.img

iput-container:
	rm -f singularity/$(APP).img.xz
	xz --compress --force --keep singularity/$(APP).img
	iput -fKP singularity/$(APP).img.xz

iget-container:
	iget -fKP $(APP).img.xz
	xz --decompress --force --keep $(APP).img.xz
	mv $(APP).img singularity/
	mv $(APP).img.xz stampede/

lytic-rsync-dry-run:
	rsync -n -arvzP --delete --exclude-from=rsync.exclude -e "ssh -A -t hpc.arizona.edu ssh -A -t lytic" ./ :project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering

lytic-rsync:
	rsync -arvzP --delete --exclude-from=rsync.exclude -e "ssh -A -t hpc.arizona.edu ssh -A -t lytic" ./ :project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering
