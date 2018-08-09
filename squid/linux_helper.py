"""
The Linux Helper module contains functionality to aid linux users to automate some tasks.

- :func:`clean_up_folder`

------------

"""

def clean_up_folder(path, files_to_remove=[], remove_empty_folders=False, verbose=False):
	"""
	Automate the removal of files from a linux system.  Given a parent directory, this will
	recursively remove specified files.

	**Parameters**

		path: *str*
			Absolute path to the parent directory to be cleaned.
		files_to_remove: *list, str, optional*
			List of files to be removed.  Wildcards can be used, thus ["*.txt","*.log"] would
			delete every file with the .txt and .log extension.
		remove_empty_folders: *bool, optional*
			Whether to remove empty folders (True), or not (False).
		verbose: *bool, optional*
			Whether to output commands used (True), or not (False).

	**Returns**

		None
	"""
	if len(files_to_remove) == 0 and remove_empty_folders is False:
		raise Exception("clean_up_folders requires either a file identifier to remove files or the remove_empty_folders flag to be set to True.")
	if path[0] != "/":
		raise Exception("For safety reasons, we require a full path to be used in clean_up_folders.")

	if path.endswith("/"): path = path[:-1]

	if verbose: print("\n---------------------\nCleaning up folder %s" % path)
	if len(files_to_remove) > 0:
		# Remove all empty folders
		ids = " -o -name ".join( map(lambda x: '"' + x + '" -delete ', files_to_remove) )
		ids = "-name " + ids
		cmd_delete_files = "find %s -type f %s" % (path, ids)
		if verbose: print("Running command: %s" % cmd_delete_files)
		os.system(cmd_delete_files)

	if remove_empty_folders:
		cmd_delete_folders = "find %s -type d -empty -delete" % path
		if verbose: print("Running command: %s" % cmd_delete_folders)
		os.system(cmd_delete_folders)
	if verbose: print("Done cleaning folder %s\n---------------------\n\n" % path)
