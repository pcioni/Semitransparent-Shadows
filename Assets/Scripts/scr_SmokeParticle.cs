using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class scr_SmokeParticle : MonoBehaviour {

    public GameObject light = null;

    GameObject faceP = null;
    GameObject shadowP = null;

	// Use this for initialization
	void Start () {
        faceP = this.transform.GetChild (0).gameObject;
        shadowP = this.transform.GetChild (1).gameObject;

        faceP.GetComponent<MeshRenderer> ().shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
        shadowP.GetComponent<MeshRenderer> ().shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.ShadowsOnly;
    }
	
	// Update is called once per frame
	void Update () {
        faceP.transform.rotation = Quaternion.LookRotation ((Camera.main.transform.position - this.transform.position)*-1.0f);
        shadowP.transform.rotation = Quaternion.LookRotation ((light.transform.position - this.transform.position)*-1.0f);
	}
}
