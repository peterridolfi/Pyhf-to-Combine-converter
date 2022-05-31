{

  "channels": [
    {
      "name": "singlechannel",
      "samples": [
        {
          "name": "sig",
          "data": [
            5.0, 10.0, 15.0
          ],
          "modifiers": [
            {
              "name": "mu",
              "type": "normfactor",
              "data": None
            }
          ]
        },
        {
          "name": "bkg",
          "data": [
            50.0, 70.0, 80.0
          ],
          "modifiers": [
            {
              "name": "uncorr_bkguncrt",
              "type": "shapesys",
              "data": [
                5.0, 10.0, 10.0
              ]
            },
            {
                "name": "corr_shape",
                "type": "histosys",
                "data": {"hi_data": [20, 15, 15], "lo_data": [10, 5, 10]}
            }
          ]
        }
      ]
    }
  ],
  "observations": [
    {
      "name": "singlechannel",
      "data": [
        53.0, 65.0, 70.0
      ]
    }
  ],
  "measurements": [
    {
      "name": "Measurement",
      "config": {
        "poi": "mu",
        "parameters": []
      }
    }
  ],
  "version": "1.0.0"
}