Developer's guide
=================

Contributing
------------
See the `Contributing Guide <https://github.com/bsc-wdc/Qdislib/blob/master/CONTRIBUTING.md>`_.


Drafting new releases
---------------------

Follow these steps when drafting a new release:

1. Ensure that the master branch is passing the tests and that the
   `latest version of the documentation <https://Qdislib.bsc.es/en/latest>`_
   is properly being built.

2. Decide whether to issue a minor or a major release following this
   `guide <https://semver.org/>`_.

3. Create and switch to a new branch named ``release-X.Y``.

4. Update the release number accordingly in the
   `VERSION <https://github.com/bsc-wdc/Qdislib/blob/master/VERSION>`_
   file.

5. Update the required PyCOMPSs version in the
   `quickstart guide <https://github.com/bsc-wdc/Qdislib/blob/master/QUICKSTART.md>`_
   if necessary.

6. Update the
   `change log <https://github.com/bsc-wdc/Qdislib/blob/master/CHANGELOG.md>`_.

7. Push the release branch with the changes.

8. Merge the newly created branch to the master branch.

9. Draft a new release in
   `Github <https://github.com/bsc-wdc/Qdislib/releases>`_ using this
   `template <https://github.com/bsc-wdc/Qdislib/blob/master/.github/RELEASE_TEMPLATE.md>`_
   using tag name ``vX.Y.Z``.

10. Create and tag a docker image for the release running the following at the
    repo's root:

    - Create the image:

      .. code:: bash

       docker build -t bscwdc/Qdislib:vX.Y.Z .

       # Create also new 'latest' tag using newly created image
       docker tag bscwdc/Qdislib:vX.Y.Z bscwdc/Qdislib:latest

    - Log in and push it to dockerhub

      .. code:: bash

       docker login -u DOCKERHUB_USER -p DOCKERHUB_PASSWORD
       docker push bscwdc/Qdislib:vX.Y.Z
       docker push bscwdc/Qdislib:latest

11. Create a pip package and upload it to PyPi:

    - Ensure that you have the latest version of ``setuptools``,
      ``wheel``, and ``twine`` installed:

      .. code:: bash

        pip3 install --upgrade setuptools wheel twine

    - Create and upload the pip package:

      .. code:: bash

       ./build.sh
       python3 -m twine upload dist/Qdislib-X.Y.Z*
