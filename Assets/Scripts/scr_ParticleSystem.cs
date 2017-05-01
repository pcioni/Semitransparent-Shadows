using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Thinksquirrel.Fluvio;

public class scr_ParticleSystem : MonoBehaviour {

    public GameObject baseParticle = null;
    FluidParticleSystem parSys = null;

    GameObject[] smokeParticles;

    int numParticles = 200;

	// Use this for initialization
	void Start () {
        parSys = this.GetComponent<FluidParticleSystem> ();
        numParticles = parSys.GetParticleSystem().main.maxParticles;

        if (baseParticle != null) {
            smokeParticles = new GameObject[numParticles];
            smokeParticles [0] = baseParticle;
            for (int i = 1; i < numParticles; i++) {
                smokeParticles [i] = GameObject.Instantiate (baseParticle);
                smokeParticles [i].SetActive (false);
            }
            smokeParticles [0].SetActive (false); 
        }

	}
	
	// Update is called once per frame
	void Update () {
        AssignParticles ();
	}

    void AssignParticles() {
        int num = parSys.GetParticleCount ();
        ParticleSystem.Particle[] particles = new ParticleSystem.Particle[num];
        parSys.GetParticleSystem ().GetParticles (particles);
        for (int i = 0; i < numParticles; i++) {
            if (i < num) {
                smokeParticles [i].SetActive (true);
                smokeParticles [i].transform.position = particles [i].position;
            }
            else {
                smokeParticles [i].SetActive (false);
            }
        }
    }

    void checkCount() {
        Debug.Log ("Number of particles: " + parSys.GetParticleSystem ().particleCount);
    }
}
