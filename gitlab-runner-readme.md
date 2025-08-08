## Starting a Gitlab Runner

First create some docker volumes for the runner to store some stuff

```
docker volume create gitlabrunner_config
docker volume create gitlabrunner_cache
```

Register the Runner with our repo:

```
docker run --rm -t -i -v gitlabrunner_config:/etc/gitlab-runner gitlab/gitlab-runner register
```

You will get asked a series of questions in order to register it.

The following is a sample of this Q&A session with provided answers marked as:
`>>like_this <<`

```
➜ docker run --rm -t -i -v gitlabrunner_config:/etc/gitlab-runner gitlab/gitlab-runner register
Runtime platform                                    arch=amd64 os=linux pid=6 revision=a8a019e0 version=12.3.0
Running in system-mode.

Please enter the gitlab-ci coordinator URL (e.g. https://gitlab.com/):
>> https://gitlab.com <<
Please enter the gitlab-ci token for this runner:
>> wq2xkxTFPESXZwow5iHV <<
Please enter the gitlab-ci description for this runner:
[4b387a33a8a9]: >> Stratos Serenity laptop <<
Please enter the gitlab-ci tags for this runner (comma separated):
>> bp-docs <<
Registering runner... succeeded                     runner=wq2xkxTF
Please enter the executor: docker-ssh, parallels, shell, ssh, docker-ssh+machine, custom, docker, virtualbox, docker+machine, kubernetes:
>> docker <<
Please enter the default Docker image (e.g. ruby:2.6):
>> alpine:latest <<
Runner registered successfully. Feel free to start it, but if it's running already the config should be automatically reloaded!
```

Keep the same answers, and maybe change the "gitlab-ci description" so we can
differentiate the runners (if need be)

Finally, start the runner with:

```
docker run --rm -t -i -v gitlabrunner_cache:/cache -v gitlabrunner_config:/etc/gitlab-runner -v /var/run/docker.sock:/var/run/docker.sock:ro gitlab/gitlab-runner
```

It will output some log entries and be ready in a sec.

You should be able to see that the runner has properly registered with our repo by
visiting: https://gitlab.com/BeyondPlanck/papers/-/settings/ci_cd, in the Runners tab (hopefully access rights allows you to do so)

To stop it, press CTLR+C and wait a few secs to wind down.

If you also have `docker-compose` installed all these start/stop commands can be much
easier with:

```
docker-compose up
```

and CTRL+C to stop it.

You can also run it as a daemon, in the background, with:

```
docker-compose up -d
```

check it's status (if it is running or not) with:

```
docker-compose ps
```

and finally stop it with:

```
docker-compose stop
```

Thanks!
